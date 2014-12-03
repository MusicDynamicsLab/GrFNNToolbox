%% canonMake
%
%  Creates "canonical" rhythm matricies, using the MIDI Toolbox note
%  representation. Rhtyhms can optionally be output to MIDI files using the
%  MIDI Toolbox 'writemidi' routine. The output can be any of all possible
%  permutations of the rhytym, in which one or more individual notes are iteratively
%  left out.
%
%  notes = canonMake(numer, denom, mType, mCode, varargin)
%
%  INPUT
%  numer - numerator of desired time signature. For complex time sigs (see
%           below), this must equal the sum of the numerator divisions.
%  denom - denominator of desired time signature
%
%  mType & mCode - metric type and code
%    Lederahl & Jackendoff metric structure specification.
%           See below for algorithm details.
%       mType = 'LJ' or 'lj'
%       mCode = vector of 1 or more ratios between metric levels, starting
%           from the lowest (shortest) levels.
%           e.g., for compound 6/8, you'd pass "[3 2]".
%           It's ok to pass irrationals for polyrhythms (e.g. 1.33), because the
%           product of it and other ratios is rounded.
%   Complex Meter
%       mType = 'CX' or 'cx'
%       mCode = vector of 1 or more numerator divisions. The sum must equal
%           the numerator of the time signature. This represent a complex
%           time signature. The first note gets a strong stress, then the begin 
%           of each division gets a weak stress, and the remainging notes are unstressed.
%           This can also be used to create regular time signatures with
%           irregular stress patterns.
%     
%  Option parameters - param/value pairs:
%
%  'outfile' -   optional filename for output. Uses MIDI Toolbox 'writemidi'.
%                   The file is written using 'numer', 'denom', 'bpm' and 'ticks' params.
%  'perm' -      [default 1] Single permutation # to create. They're made by iteratively
%                           leaving out one or more notes. The numbering comes from the
%                           algorithm used to create the permutations.
%                           Pass 1 for full rhythm. If out-of-range, returns full rhythm.
%  'permCount' - [default 0] Flag. Pass 1 to simply return the # of
%                           permutations possible for the passed
%                           parameters. This takes the 'permIncludeFirst' & 'permMin' 
%                           flags into account.
%                           The permutation matrix is displayed.
%  'permAll' -   [default 0] Flag. Pass 1 to have all permutations created. 
%                       A cell matrix is returned containing all permutations. 
%                       If an outfile is specified, each is also written to an outfile, 
%                       and each permutation # is appended to outfile name.
%  'permMin' -   [default 2] - Minimum # of notes to have in a permutation. 
%                           This is to avoid permutations simply with a
%                           note on first beat.
%                           With value 2 and permIncludeFirst = 0, the 
%                           most sparse permutations will 
%                           have a note on 1st beat, and then one other position.
%  'permIncludeFirst' - [default 0] Flag. Pass 1 to have the permutations include
%                           *leaving out* the first note, i.e. the downbeat.
%  'permV'   -   Permutation Vector. Pass a row vector of 1's and 0's to
%                   specify a permutation mask to use instead of choosing
%                   by 'perm' param. Must be of proper length, i.e. as long
%                   as single repitition of the full rhythm.
%  'stress' -    [default 'auto'] stress pattern
%     one of:
%     'none' -  all notes unstressed
%     'first' - stress first beat/note only
%     'auto'    [default] - stress strong and weak notes, according to
%                           metric structure.
%  Stress velocities:
%    'strong' -  [default 127] MIDI velocity of strong stress, for first beat
%    'weak' -    [default 96] MIDI velocity of weak stress, used in stress type 'auto'.
%    'unstressed' - [default 64] MIDI velocity of unstressed beats.
%
%  'pitch' -     [default 48] MIDI pitch value for notes
%  'reps'   -    [default 1] Total # of repetitions of the rhythm. >= 1.
%  'endSilence'  [default 0] Flag. Pass 1 to have N measures of silence after
%                           the N repitions of the rhythm, with a single note 
%                           at begin of measure 2N+1 to signal the end.
%  'bpm' -       [default 120] - bpm, used to calculate event realtime second values, 
%                           and passed to optional file output.
%                           *NOTE* This always specifies
%                           quarter-notes/minute, regardless of time signature.
%  'ticks' -     [default 120] - ticks per beat for optional file output.
%                               12/15/08 Not recommended to use something other than
%                               120. Found a bug recently in MIDI Toolbox writemidi
%                               that gave some buggy output if ticks-per-beat is not
%                               120. Not sure it would effect these small files here,
%                               but better to be safe.
%
%
%  OUTPUT
%   Returns
%    If 'permCount' option is set to 1:
%       - displays the permutation matrix
%       - returns the # of possible permutations for the passed settings.
%
%
%  Otherwise:
%    notes - matrix of notes in 'nmat' format used by MIDI Toolbox.
%           Each row represents one note, with these values by column:
%               onset(beats), duration(beats), MIDI Channel, MIDI Pitch,
%                   MIDI velocity, onset(seconds), duration(seconds)
%           If 'permAll' is passed, output is a cell matrix with each permutation.
%
%           NOTE that 'beats' actually means quarter notes. So that
%           regardless of time signature, the beats are in terms of quarter
%           notes. This is how MIDI Toolbox does it. Also, the onset and
%           duration in terms of seconds are calculated using the bpm in
%           terms of quarter-notes per minute, regardless of time
%           signature.
%
%  DETAILS
%    Lederal & Jackendoff metric structures
%       The vector of ratios is used to construct the rhythm in conjunction
%       with the time signature. The time signature determines the duration
%       in quarter notes of the full rhythm (a single measure). Each metric
%       level contributes N equally-spaced notes within the rhythm, where N is # of
%       notes on that level, deduced from the vector of ratios. These notes
%       get a weak stress, except for the lowest (shortest) level, which is
%       unstressed. This allows for defining polyrhythms, such as [1.5 2],
%       which will yield 4 notes, with a strong stress on 1 and weak on the 
%       3rd note ('2-and' if 3/4 time signature).
%       The first note always gets a strong stress.
%
%       Although some polyrhythms will seem to be duplicated by permutations
%       of straight rhythms, the stress patterns of polyrhythms are
%       different, and thus can be useful.

%%

function notes = canonMake(numer, denom, mType, mCode, varargin)


%% Params - check, get and set defaults.
ps = inputParser;
ps.FunctionName = 'canonMake';
ps.addRequired('numer', @(x)x>1 && mod(x,1)==0);
ps.addRequired('denom', @(x)numel(find([2 4 8 16 32] == x))>0);
ps.addRequired('mType', @ischar );
ps.addRequired('mCode', @isnumeric);
ps.addParamValue('outfile','', @ischar);
ps.addParamValue('perm', 1, @(x)x>0&&isnumeric(x));
ps.addParamValue('permCount', 0, @(x)x>=0&&isnumeric(x));
ps.addParamValue('permAll',0,@(x)x>=0&&isnumeric(x));
ps.addParamValue('permMin',2,@(x)x>0&&isnumeric(x));
ps.addParamValue('permIncludeFirst',0,@(x)x>=0&&isnumeric(x));
ps.addParamValue('permV','',@isnumeric);
ps.addParamValue('stress', 'auto', @ischar);
ps.addParamValue('strong', 127, @(x) isnumeric(x) && x>=0 && x<128);
ps.addParamValue('weak', 96, @(x) isnumeric(x) && x>=0 && x<128);
ps.addParamValue('unstressed', 64, @(x) isnumeric(x) && x>=0 && x<128);
ps.addParamValue('pitch', 48, @(x) isnumeric(x) && x>=1 && x<128);
ps.addParamValue('reps', 1, @(x)x>0&&isnumeric(x));
ps.addParamValue('endSilence', 0, @(x)x>=0&&isnumeric(x));
ps.addParamValue('bpm', 120, @isnumeric);
ps.addParamValue('ticks',120,@(x)x>=120&&isnumeric(x));

ps.parse(numer, denom, mType, mCode, varargin{:});
p = ps.Results;

%% check ticks value
if( p.ticks ~= 120)
    warning('\ncanonMake: ticks value is not 120. Read warnings in code comments. Proceeding.\n'); %#ok<WNTAG>
end

%% Modify stress settings
if strcmpi(p.stress,'none')
    %all notes unstressed
    p.strong = p.unstressed;
    p.weak = p.unstressed;
elseif strcmpi(p.stress,'first')
    %only first beat is stressed, so set weak stress value to unstressed
    p.weak = p.unstressed;
elseif strcmpi(p.stress,'auto') == 0
    error('unrecognized stress param value');
    return;
end

%% Ticks per quarter note, used for internal calcs.
%  Realworld timing is done using bpm = quarter note / minute
tpq = 120; %NOTE: see comments above
%ticks per measure (single loop of rhythm). denom is {2,4,8,16,32}
tpm = p.numer * tpq * 4 / p.denom;

%output note matrix
notes = [];


%% Complex meter
%
if strcmpi(p.mType,'CX')
    %The numerator of the input time signature must match the 
    % sum of the divisions in the complex time signature vector.
    if sum(p.mCode) ~= p.numer
        error('canonMake: sum of complex time signature divisions must equal numerator.');
        return;
    end
    dur = tpq * 4 / p.denom;
    inc = 0;
    for d = 1:length(p.mCode)
       for n = 1:p.mCode(d)
           vel = p.unstressed;
           if n == 1 %begin of each division is weakly stressed
               vel = p.weak;
           end
           notes = [notes; [dur*inc 0 1 p.pitch vel 0 0] ];
           inc = inc+1;
       end
    end
%%%  

%% Lederal & Jackendoff
%
elseif strcmpi(p.mType,'LJ')
    %Go through each LJ level-ratio in the mCode vector and add notes at the
    % appropriate note duration for that level. This will handle polyrhythm
    % correctly.
    numNotes = 1;
    %Loop from highest (longest) level. The vector goes from lowest to
    %highest.
    for l = length(p.mCode):-1:1
       numNotes = round(numNotes .* mCode(l));
            %call round() cuz for polyrhtyhms, can have ratios of 1.33
       dur = tpm / numNotes;
       vel = p.weak;
       if l == 1
           %The lowest(shortest) level is not stressed (except first beat, handled below)
           vel = p.unstressed;
       end
       for nn=1:numNotes
          notes = [notes; [dur*(nn-1) 0 1 p.pitch vel 0 0] ];
       end
    end
    %Sort and remove duplicates.
    [b,i,j]=unique(notes(:,1),'first'); 
    %i holds array of first occurences of unique rows. Preserves
    % velocities of higher levels.
    notes = notes(i,:);
else
    error('canonMake: unrecognized metric type');
    return;
end    

%% Check if want back is the # of permutations.
if p.permCount == 1   
   prms = calcPerms(size(notes,1), p.permIncludeFirst, p.permMin);
   nums = 1:1:size(prms,1);
   disp('     #     Permutation Mask');
   disp([nums' prms]);
   notes = size(prms,1);
   return;
end

%% Caclulate durations
for l=1:size(notes,1)
    if l == size(notes,1)
       %last note. make it last till end of measure
       notes(l,2) = tpm - notes(l,1); %#ok<AGROW>
    else
       %otherwise, it lasts til next note. Handles polyrhythms.
       notes(l,2) = notes(l+1,1) - notes(l,1); %#ok<AGROW>
    end
end

%Full stress on first note
notes(1,5) = p.strong; %is set to unstressed value when 'none' passed for stress type

%% Check Permutations
%  Check for a permutation vector to use directly
if ~isempty(p.permV)
    %Check length
    if( length(p.permV) ~= size(notes,1) )
        error('canonMake: permV vector must be length of rhythm: %d',size(notes,1));
    end
    notes = notes(find(p.permV == 1), :); %#ok<FNDSB>
    
% Determine which permutation to get, or all.
elseif p.permAll == 1 || p.perm > 1
   %Get array of all possible permutations, as rows of 1's && 0's, representing which
   % notes should be included in each perm.
   prms = calcPerms(size(notes,1), p.permIncludeFirst, p.permMin);
   if p.permAll == 1
      %return a cell array of each permutation.
      notesCell = {};
      for n=1:size(prms,1)
         notesCell{n} = notes(find(prms(n,:) == 1), :); %#ok<FNDSB,AGROW>
      end
      notes = notesCell;
   else
      %Return a single permutation
      %If the perm # is out of range, just return full rhythm.
      %We already tested that it's > 0 during init.
      if p.perm <= size(prms,1)
        notes = notes(find(prms(p.perm,:) == 1), :); %#ok<FNDSB>
      else
        warning('perm parameter out-of-range. Returning full rhythm');  %#ok<WNTAG>
      end
   end
end

%% Finish up
num = 1;
if p.permAll == 1
    num = length(notesCell);
end

for nn = 1:num
    if p.permAll == 1
        notes = notesCell{nn};
    end
    
    %Do repetitions
    if p.reps > 1
       orig = notes;
       for r = 1:p.reps-1
          next = orig;
          next(:,1) = next(:,1) + r*tpm;
          notes = [notes; next];  %#ok<AGROW>
       end
    end

    %Add end-note for silence
    if p.endSilence > 0
        notes = [notes; [tpm*(2*p.reps) tpq 1 p.pitch p.unstressed 0 0]]; %#ok<AGROW>
    end

    %Convert from ticks to beat (quarter-note, actually) numbering
    notes(:,[1 2]) = notes(:,[1 2]) / tpq;

    %Calc the realtime onsets and durations. 
    notes(:,6) = notes(:,1) * 60 / p.bpm;
    notes(:,7) = notes(:,2) * 60 / p.bpm;

    if p.permAll == 1
        notesCell{nn} = notes; %#ok<AGROW>
        notes = notesCell;
    end
end

%Output file
if numel(p.outfile) > 0
   if p.permAll == 1
       %Write out each permutation to a separate file.
       [pathstr, name, ext, versn] = fileparts(p.outfile); %#ok<NASGU>
       for n=1:length(notesCell)
            file = fullfile(pathstr,[name sprintf('%05d',n) ext]);
            writemidi(notesCell{n}, file, p.ticks, p.bpm, p.numer, p.denom);
       end
   else
       writemidi(notes, p.outfile, p.ticks, p.bpm, p.numer, p.denom);
   end
end

end

%% Utility Functions

%%%%%%%%%%%%%%%%%%%%%%%%%
%Permutations. 
%Creates an array of vectors, each row corresponding to a possible
% permutation of 1 and 0, with length N.
%Permutations are based on whether or not a given note from the full
% rhythm is included in the output, so it's a binary permutation.
%Do a find(res(n)) to find which notes should be included in 
% permutation n.
%This is much more memory efficient than using builtin "perms" func,
% cuz that calcs all possible permutations by position which results
% in massive amount of permutations when you get up around N = 10, 11.
%
function res = calcPerms(N, includeFirstBeat, permMin)
    res = zeros(2^N, N);
    res = bifur(res);
    %Check if should remove perms that start with 0
    if(includeFirstBeat == 0)
        res = res(1:end/2, :);
    end
    %Remove perms that don't have enough 1's (i.e., enough notes).
    for n=1:size(res,1)
       keep(n) = numel(find(res(n,:) > 0)) >= permMin;
    end
    res = res(find(keep > 0),:);
end

function res = bifur(res)
    %"bifurcate" the passed array, setting '1' to first column of top half,
    % and '0' to first column of bottom half. Then call again sub-array, to
    % recursively create the whole permutations array.
    res(1:end/2,1) = 1;
    res(end/2+1:end,1) = 0;
    if size(res,2) > 1
        res(1:end/2, 2:end) = bifur(res(1:end/2,2:end));
        res(end/2+1:end, 2:end) = bifur(res(end/2+1:end, 2:end));
    end
end

%% modifications
%  ** MODIFICATIONS **
%  01/07/08  MGS Created.
%  12/15/08  MGS Changed default 'ticks' value to 120. Added warning if
%                not 120. See comments above under 'ticks' param.
%  01/12/09  MGS Changed 'endSilence' default to 0. Intending for main use
%                to come from making sub-rhythms for subsequent concatenation.
%                Changed 'reps' default to 1. Like above for 'endSilence'.
%                Improved permutation display when set 'permCount';
%                Adding 'permV' option to pass perumation mask vector
%                directly.