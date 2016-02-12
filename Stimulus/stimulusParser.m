%% stimulusParser
% Error checker and interpolator for all user inputs to stimulusMake. Not a standalone function.
%

%%
function s = stimulusParser(s, varargin)

if strcmpi(s.type,'fcn')
    
    %% Time courses and frequencies
    
    if size(s.ts,1) > 1 && size(s.fc,1) == 1
        s.fc = repmat(s.fc,size(s.ts,1),1);
        %     warning('Frequency or complex copied to all time courses since only one was specified')
    elseif size(s.fc,1) > size(s.ts,1)
        error('More frequency time courses specified than time courses')
    elseif size(s.fc,1) < size(s.ts,1)
        error('Less frequency time courses specified than time courses')
    end
    if ~iscell(s.fc), s.fc = num2cell(s.fc);end
    
    %% Initial phases
    
    s.Th = zeros(size(s.fc));
    
    if length(varargin) > 5 && ~ischar(varargin{6})
        s.Th = varargin{6};     %initial phases of components, written here if spec'd in varargin
    end
    
    if numel(s.Th) == 1 && numel(s.fc) ~= 1            %commence error checking for sizes, # of attribute values, etc.
        s.Th = repmat(s.Th,size(s.fc,1),size(s.fc,2));
        %     warning('Initial phase copied to all components and/or time courses since only one was specified')
    elseif size(s.Th,1) == 1 && size(s.fc,1) > 1 && size(s.Th,2) == size(s.fc,2)
        s.Th = repmat(s.Th,size(s.fc,1),1);
        %     warning('Initial phase(s) for each component copied to all time courses since only one was specified')
    elseif size(s.Th,2) == 1 && size(s.fc,2) > 1 && size(s.Th,1) == size(s.fc,1)
        s.Th = repmat(s.Th,1,size(s.fc,2));
        %     warning('Initial phase for each time course copied to all components since only one was specified for each')
    elseif size(s.Th,1) ~= size(s.fc,1)
        error('Must specify initial phase for every component of every timecourse, or input just a scalar or vector that will be copied')
    elseif size(s.Th,2) ~= size(s.fc,2)
        error('Must specify initial phase for every component of every timecourse, or input just a scalar or vector that will be copied')
    end
    
    %% Carrier waveforms
    
    if numel(s.carrier) == 1 && numel(s.fc) ~= 1
        s.carrier = repmat(s.carrier,size(s.fc,1),size(s.fc,2));
        %     warning('Carrier waveform copied to all components and/or time courses since only one was specified')
    elseif size(s.carrier,1) == 1 && size(s.fc,1) > 1 && size(s.carrier,2) == size(s.fc,2)
        s.carrier = repmat(s.carrier,size(s.fc,1),1);
        %     warning('Carrier waveforms for components copied to all time courses since only one was specified')
    elseif size(s.carrier,2) == 1 && size(s.fc,2) > 1 && size(s.carrier,1) == size(s.fc,1)
        s.carrier = repmat(s.carrier,1,size(s.fc,2));
        %     warning('Carrier waveform for each time course copied to all components since only one was specified for each')
    elseif size(s.carrier,1) ~= size(s.fc,1)
        error('Must specify carrier waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
    elseif size(s.carrier,2) ~= size(s.fc,2)
        error('Must specify carrier waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
    end
    
    %% Amplitudes
    
    if numel(s.ac) == 1 && numel(s.fc) ~= 1
        s.ac = repmat(s.ac,size(s.fc,1),size(s.fc,2));
        %     warning('Amplitude copied to all components and/or time courses since only one was specified')
    elseif size(s.ac,1) == 1 && size(s.fc,1) > 1 && size(s.ac,2) == size(s.fc,2)
        s.ac = repmat(s.ac,size(s.fc,1),1);
        %     warning('Amplitudes for components copied to all time courses since only one was specified')
    elseif size(s.ac,2) == 1 && size(s.fc,2) > 1 && size(s.ac,1) == size(s.fc,1)
        s.ac = repmat(s.ac,1,size(s.fc,2));
        %     warning('Amplitude for each time course copied to all components since only one was specified for each')
    elseif size(s.ac,1) ~= size(s.fc,1)
        error('Must specify amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
    elseif size(s.ac,2) ~= size(s.fc,2)
        error('Must specify amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
    end
    if ~iscell(s.ac), s.ac = num2cell(s.ac);end
    
    %% Check if correct number of values for each additional attribute
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i},'ramp')
            if length(varargin) > i + 1 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2})
                s.sc = varargin{i+1}; %ramp time in ms of beginning and end of stim
                s.sp = varargin{i+2}; %ramp strength exponent, ie, larger than 1, sudden onset, smaller than one, gradual
                if length(varargin) > i + 2 && ~ischar(varargin{i+3})
                    error('Ramp takes two variables: Length in ms and exponent of curve')
                end
            else
                error('Ramp takes two variables: Length in ms and exponent of curve')
            end
        end
        if strcmpi(varargin{i},'fm')
            if length(varargin) > i + 2 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3})
                s.FM = varargin{i+1};      % FM modulator waveforms
                s.fFM = varargin{i+2};     % FM modulation frequencies
                s.aFM = varargin{i+3};     % FM modulation amplitudes
                if length(varargin) > i + 3 && ~ischar(varargin{i+4})
                    error('FM takes three variables: Waveform, frequency and amplitude (depth)')
                end
            else
                error('FM takes three variables: Waveform, frequency and amplitude (depth)')
            end
        end
        if strcmpi(varargin{i},'am')
            if length(varargin) > i + 3 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2})...
                    && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4})
                s.AM = varargin{i+1};      % AM modulator waveforms
                s.fAM = varargin{i+2};     % AM modulation frequencies
                s.aAM = varargin{i+3};     % AM modulation depths
                s.CfreqsAM = varargin{i+4};% ampl. of carrier freqs remaining in stim (0 to 1)
                if length(varargin) > i + 4 && ~ischar(varargin{i+5})
                    error('AM takes four variables: Waveform, frequency, amplitude (depth) and ampl. of carrier in final wave')
                end
            else
                error('AM takes four variables: Waveform, frequency, amplitude (depth) and ampl. of carrier in final wave')
            end
        end
        if strcmpi(varargin{i},'mask')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.mask = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Noise (for masking stim time course(s)) takes one variable: RMS signal-to-noise ratio(s) (in dB)')
                end
            else
                error('Noise (for masking stim time course(s)) takes one variable: RMS signal-to-noise ratio(s) (in dB)')
            end
        end
        if strcmpi(varargin{i},'iter')
            if length(varargin) > i + 1 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2})
                s.iter = varargin{i+1};
                s.Niter = varargin{i+2};
                if length(varargin) > i + 2 && ~ischar(varargin{i+3})
                    error('Iter takes two variables: Length in seconds of delay and number of iterations')
                end
            else
                error('Iter takes two variables: Length in seconds of delay and number of iterations')
            end
        end
        if strcmpi(varargin{i},'filtstim')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtstim = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
                end
            else
                error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
            if size(s.filtstim,2) ~= 2
                error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
        end
        if strcmpi(varargin{i},'filtmask')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtmask = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
                end
            else
                error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
            if size(s.filtmask,2) ~= 2
                error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
        end
        if strcmpi(varargin{i},'maskall')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.maskall = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
                end
            else
                error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
            end
            if numel(s.maskall) ~= 1
                error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
            end
        end
        if strcmpi(varargin{i},'filtmaskall')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtmaskall = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
                end
            else
                error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
            end
            if numel(s.filtmaskall) ~= 2
                error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
            end
        end
    end
    
    %% Ramp
    
    if numel(s.sc) == 1 && size(s.fc,1) ~= 1
        s.sc = repmat(s.sc,size(s.fc,1),1);
        %     warning('Ramp time copied to all time courses since only one was specified')
    elseif numel(s.sc) > 1 && size(s.fc,1) ~= numel(s.sc)
        error('Must specify ramp time for every timecourse, or input just a scalar that will be copied')
    end
    
    if numel(s.sp) == 1 && size(s.fc,1) ~= 1
        s.sp = repmat(s.sp,size(s.fc,1),1);
        %     warning('Ramp exponent copied to all time courses since only one was specified')
    elseif numel(s.sp) > 1 && size(s.fc,1) ~= numel(s.sp)
        error('Must specify ramp exponent for every timecourse, or input just a scalar that will be copied')
    end
    
    %% Frequency modulator
    
    if isfield(s, 'FM')
        if numel(s.fFM) == 1 && numel(s.fc) ~= 1
            s.fFM = repmat(s.fFM,size(s.fc,1),size(s.fc,2));
            %         warning('FM freq copied to all components and/or time courses since only one was specified')
        elseif size(s.fFM,1) == 1 && size(s.fc,1) > 1 && size(s.fFM,2) == size(s.fc,2)
            s.fFM = repmat(s.fFM,size(s.fc,1),1);
            %         warning('FM freqs for each component copied to all time courses since only one was specified')
        elseif size(s.fFM,2) == 1 && size(s.fc,2) > 1 && size(s.fFM,1) == size(s.fc,1)
            s.fFM = repmat(s.fFM,1,size(s.fc,2));
            %         warning('FM freq for each time course copied to all components since only one was specified for each')
        elseif size(s.fFM,1) ~= size(s.fc,1)
            error('Must specify FM freq for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.fFM,2) ~= size(s.fc,2)
            error('Must specify FM freq for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
        if ~iscell(s.fFM), s.fFM = num2cell(s.fFM);end
        
        if numel(s.FM) == 1 && numel(s.fc) ~= 1
            s.FM = repmat(s.FM,size(s.fc,1),size(s.fc,2));
            %         warning('FM waveform copied to all components and/or time courses since only one was specified')
        elseif size(s.FM,1) == 1 && size(s.fc,1) > 1 && size(s.FM,2) == size(s.fc,2)
            s.FM = repmat(s.FM,size(s.fc,1),1);
            %         warning('FM waveforms for each component copied to all time courses since only one was specified')
        elseif size(s.FM,2) == 1 && size(s.fc,2) > 1 && size(s.FM,1) == size(s.fc,1)
            s.FM = repmat(s.FM,1,size(s.fc,2));
            %         warning('FM waveform for each time course copied to all components since only one was specified for each')
        elseif size(s.FM,1) ~= size(s.fc,1)
            error('Must specify FM waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
        elseif size(s.FM,2) ~= size(s.fc,2)
            error('Must specify FM waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
        end
        
        if numel(s.aFM) == 1 && numel(s.fc) ~= 1
            s.aFM = repmat(s.aFM,size(s.fc,1),size(s.fc,2));
            %         warning('FM amplitude copied to all components and/or time courses since only one was specified')
        elseif size(s.aFM,1) == 1 && size(s.fc,1) > 1 && size(s.aFM,2) == size(s.fc,2)
            s.aFM = repmat(s.aFM,size(s.fc,1),1);
            %         warning('FM amplitudes for each component copied to all time courses since only one was specified')
        elseif size(s.aFM,2) == 1 && size(s.fc,2) > 1 && size(s.aFM,1) == size(s.fc,1)
            s.aFM = repmat(s.aFM,1,size(s.fc,2));
            %         warning('FM amplitude for each time course copied to all components since only one was specified for each')
        elseif size(s.aFM,1) ~= size(s.fc,1)
            error('Must specify FM amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.aFM,2) ~= size(s.fc,2)
            error('Must specify FM amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
        if ~iscell(s.aFM), s.aFM = num2cell(s.aFM);end
    end
    
    %% Amplitude modulator
    
    if isfield(s, 'AM')
        if numel(s.fAM) == 1 && numel(s.fc) ~= 1
            s.fAM = repmat(s.fAM,size(s.fc,1),size(s.fc,2));
            %         warning('AM freq copied to all components and/or time courses since only one was specified')
        elseif size(s.fAM,1) == 1 && size(s.fc,1) > 1 && size(s.fc,2) == size(s.fAM,2)
            s.fAM = repmat(s.fAM,size(s.fc,1),1);
            %         warning('AM freqs for each component copied to all time courses since only one was specified')
        elseif size(s.fAM,2) == 1 && size(s.fc,2) > 1 && size(s.fc,1) == size(s.fAM,1)
            s.fAM = repmat(s.fAM,1,size(s.fc,2));
            %         warning('AM freq for each time course copied to all components since only one was specified for each')
        elseif size(s.fAM,1) ~= size(s.fc,1)
            error('Must specify AM freq for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.fAM,2) ~= size(s.fc,2)
            error('Must specify AM freq for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
        if ~iscell(s.fAM), s.fAM = num2cell(s.fAM);end
        
        if numel(s.AM) == 1 && numel(s.fc) ~= 1
            s.AM = repmat(s.AM,size(s.fc,1),size(s.fc,2));
            %         warning('AM waveform copied to all components and/or time courses since only one was specified')
        elseif size(s.AM,1) == 1 && size(s.fc,1) > 1 && size(s.AM,2) == size(s.fc,2)
            s.AM = repmat(s.AM,size(s.fc,1),1);
            %         warning('AM waveforms for each component copied to all time courses since only one was specified')
        elseif size(s.AM,2) == 1 && size(s.fc,2) > 1 && size(s.AM,1) == size(s.fc,1)
            s.AM = repmat(s.AM,1,size(s.fc,2));
            %         warning('AM waveform for each time course copied to all components since only one was specified for each')
        elseif size(s.AM,1) ~= size(s.fc,1)
            error('Must specify AM waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
        elseif size(s.AM,2) ~= size(s.fc,2)
            error('Must specify AM waveform for every component of every timecourse, or input just a single one or a vector that will be copied')
        end
        
        if numel(s.aAM) == 1 && numel(s.fc) ~= 1
            s.aAM = repmat(s.aAM,size(s.fc,1),size(s.fc,2));
            %         warning('AM amplitude copied to all components and/or time courses since only one was specified')
        elseif size(s.aAM,1) == 1 && size(s.fc,1) > 1 && size(s.aAM,2) == size(s.fc,2)
            s.aAM = repmat(s.aAM,size(s.fc,1),1);
            %         warning('AM amplitudes for each component copied to all time courses since only one was specified')
        elseif size(s.aAM,2) == 1 && size(s.fc,2) > 1 && size(s.aAM,1) == size(s.fc,1)
            s.aAM = repmat(s.aAM,1,size(s.fc,2));
            %         warning('AM amplitude for each time course copied to all components since only one was specified for each')
        elseif size(s.aAM,1) ~= size(s.fc,1)
            error('Must specify AM amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.aAM,2) ~= size(s.fc,2)
            error('Must specify AM amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
        if ~iscell(s.aAM), s.aAM = num2cell(s.aAM);end
        
        if numel(s.CfreqsAM) == 1 && numel(s.fc) ~= 1
            s.CfreqsAM = repmat(s.CfreqsAM,size(s.fc,1),size(s.fc,2));
            %         warning('AM carrier amplitude copied to all components and/or time courses since only one was specified')
        elseif size(s.CfreqsAM,1) == 1 && size(s.fc,1) > 1 && size(s.CfreqsAM,2) == size(s.fc,2)
            s.CfreqsAM = repmat(s.CfreqsAM,size(s.fc,1),1);
            %         warning('AM carrier amplitudes for each component copied to all time courses since only one was specified')
        elseif size(s.CfreqsAM,2) == 1 && size(s.fc,2) > 1 && size(s.CfreqsAM,1) == size(s.fc,1)
            s.CfreqsAM = repmat(s.CfreqsAM,1,size(s.fc,2));
            %         warning('AM carrier amplitude for each time course copied to all components since only one was specified for each')
        elseif size(s.CfreqsAM,1) ~= size(s.fc,1)
            error('Must specify AM carrier amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.CfreqsAM,2) ~= size(s.fc,2)
            error('Must specify AM carrier amplitude for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
    end
    
    %% Mask noise
    
    if isfield(s, 'mask')
        if size(s.mask,1) ~= 1 && size(s.mask,2) ~= 1
            error('Mask noise must be a scalar applying to every time course or a length(# time courses) vector with values corresponding to each time course')
        elseif size(s.mask,2) == 1 && size(s.mask,1) == size(s.fc,1) && numel(s.mask) ~= 1
            s.mask = s.mask';
        elseif size(s.mask,1) == 1 && size(s.mask,2) ~= size(s.fc,1) && numel(s.mask) ~= 1
            error('Mask noise must be a scalar applying to every time course or a length(# time courses) vector with values corresponding to each time course')
        elseif numel(s.mask) == 1
            s.mask = repmat(s.mask,1,size(s.fc,1));
        end
    end
    
    %% Filters
    
    if isfield(s, 'filtstim')
        if size(s.filtstim,1) == 1 && size(s.fc,1) ~= 1
            s.filtstim = repmat(s.filtstim,size(s.fc,1),1);
            %         warning('Stimulus filter applied to all time courses since only one was specified')
        elseif size(s.filtstim,1) ~= size(s.fc,1)
            error('Number of stimulus filters must equal number of time courses, or input just one that will be copied')
        end
    end
    
    if isfield(s, 'filtmask')
        if size(s.filtmask,1) == 1 && size(s.fc,1) ~= 1
            s.filtmask = repmat(s.filtmask,size(s.fc,1),1);
            %         warning('Mask noise filter applied to all time courses since only one was specified')
        elseif size(s.filtmask,1) ~= size(s.fc,1)
            error('Number of mask noise filters must equal number of time courses, or input just one that will be copied')
        end
    end
    
    %% Delay-and-add iteration
    
    if isfield(s, 'iter')
        if numel(s.iter) == 1 && numel(s.fc) ~= 1
            s.iter = repmat(s.iter,size(s.fc,1),size(s.fc,2));
            %         warning('Delay-Add interval copied to all components and/or time courses since only one was specified')
        elseif size(s.iter,1) == 1 && size(s.fc,1) > 1 && size(s.iter,2) == size(s.fc,2)
            s.iter = repmat(s.iter,size(s.fc,2),1);
            %         warning('Delay-Add intervals for each component copied to all time courses since only one was specified')
        elseif size(s.iter,2) == 1 && size(s.fc,2) > 1 && size(s.iter,1) == size(s.fc,1)
            s.iter = repmat(s.iter,1,size(s.fc,2));
            %         warning('Delay-Add interval for each time course copied to all components since only one was specified for each')
        elseif size(s.iter,1) ~= size(s.fc,1)
            error('Must specify Delay-Add interval for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.iter,2) ~= size(s.fc,2)
            error('Must specify Delay-Add interval for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
        
        if numel(s.Niter) == 1 && numel(s.fc) ~= 1
            s.Niter = repmat(s.Niter,size(s.fc,1),size(s.fc,2));
            %         warning('Number of delay iterations copied to all components and/or time courses since only one was specified')
        elseif size(s.Niter,1) == 1 && size(s.fc,1) > 1 && size(s.Niter,2) == size(s.fc,2)
            s.Niter = repmat(s.Niter,size(s.fc,1),1);
            %         warning('Numbers of delay iterations for each component copied to all time courses since only one was specified')
        elseif size(s.Niter,2) == 1 && size(s.fc,2) > 1 && size(s.Niter,1) == size(s.fc,1)
            s.Niter = repmat(s.Niter,1,size(s.fc,2));
            %         warning('Number of delay iterations for each time course copied to all components since only one was specified for each')
        elseif size(s.Niter,1) ~= size(s.fc,1)
            error('Must specify number of delay iterations for every component of every timecourse, or input just a scalar or vector that will be copied')
        elseif size(s.Niter,2) ~= size(s.fc,2)
            error('Must specify number of delay iterations for every component of every timecourse, or input just a scalar or vector that will be copied')
        end
    end
    
elseif strcmpi(s.type,'wav')
    
    s0 = 0;
    sf = length(s.x)-1;
    s.ts = [s0 sf]/s.fs;
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i}(1:3),'ste')
            if size(s.x,2) == 1                            % Expand to stereo if mono
                s.x = [s.x s.x];
            end
        end
        if strcmpi(varargin{i}(1:3),'mon')
            if size(s.x,2) == 2
                s.x = sum(s.x,2) / size(s.x,2);            % Contract to mono if stereo
            end
        end
    end
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i},'mask')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.mask = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Noise (for masking stim time course(s)) takes one variable: RMS signal-to-noise ratio (in dB)')
                end
            else
                error('Noise (for masking stim time course(s)) takes one variable: RMS signal-to-noise ratio (in dB)')
            end
        end
        if strcmpi(varargin{i},'filtstim')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtstim = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
                end
            else
                error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
            if size(s.filtstim,2) ~= 2
                error('Filtstim takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
        end
        if strcmpi(varargin{i},'filtmask')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtmask = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
                end
            else
                error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
            if size(s.filtmask,2) ~= 2
                error('Filtmask takes one variable: An m-by-2 array, m is # of time courses, each row is vector inputs b and a to filter() function')
            end
        end
        if strcmpi(varargin{i},'maskall')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.maskall = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
                end
            else
                error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
            end
            if numel(s.maskall) ~= 1
                error('Maskall takes one scalar variable: RMS signal-to-noise ratio (in dB)')
            end
        end
        if strcmpi(varargin{i},'filtmaskall')
            if length(varargin) > i && ~ischar(varargin{i+1})
                s.filtmaskall = varargin{i+1};
                if length(varargin) > i + 1 && ~ischar(varargin{i+2})
                    error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
                end
            else
                error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
            end
            if numel(s.filtmaskall) ~= 2
                error('Filtmaskall takes one variable: A 1-by-2 array, b and a, vectors to filter() function')
            end
        end
        if strcmpi(varargin{i},'ts')
            s.newTS = varargin{i+1};
        end
        if strcmpi(varargin{i},'fs')
            s.newFS = varargin{i+1};
        end
        if strcmpi(varargin{i}(1:3),'gam')
            s.gam.minCF    = varargin{i+1};
            s.gam.maxCF    = varargin{i+2};
            s.gam.numChans = varargin{i+3};
        end
        if strcmpi(varargin{i},'ramp')
            s.sc = varargin{i+1};
            s.sp = varargin{i+2};
        end
        
    end
    
end
%% revision history
% KDL - 1/30/14 - Commenting out all warnings
