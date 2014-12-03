%% connectMake
%  C = connectMake(n1, n2, type, varargin)
% 
%  Creates a connection matrix from network n1 to network n2
%  'type' is the type of connectivity
%
%  n1                connection From here
%  n2                connection To here
%  type              'full', 'near', 'one', 'gaus'
% 
%  Optional inputs that must come in this order after type:
% 
%  amplitude         One scalar afterward, strength
%  range             Another scalar (freq ratio for log fspac, freq diff for lin)
%                    Use di2df to specify range in number of oscillators
%                    This specifies the range of frequencies type 'gaus' or
%                    'one' applies
% 
%  Other optional inputs, can come in any order:
% 
%  'complex'           No input afterwards; this will specify that
%                      connections will have both amplitude and phase
%  'modes'             Input scalar or vector afterwards, location (freq ratio) of multiple modes
%  'amps'              Input scalar or vector afterwards, relative amplitude of multiple modes
%  'ranges'            Input scalar or vector afterwards, range of multiple modes (in fraction of range)
% 
%  Example calls:
%
%   C = connectMake(n1, n2, 'gaus', .7, 1.05, 'modes', [1/2 1 2], 'amps', [.25 1 .5], 'ranges', [.3 1 .7]);
%   C = connectMake(n1, n2, 'one', .4, di2df(n1, 5), 'complex');
%   C = connectMake(n1,n2,'full',10);
%

%%
function C = connectMake(n1, n2, type, varargin)

if nargin < 3 || ~ischar(type)
    error('Must specify kernel type: full, near, one, or gaus');
end
if ~any(strcmpi(type,{'full' 'near' 'one' 'gaus'}))
    error('Must specify kernel type: full, near, one, or gaus');
end

amp  = 1;                   % Set defaults
switch n2.fspac
    case 'log'
        range = 1.02; % freq ratio
    case 'lin'
        range = .02; % freq difference
end
complexKernel = 0;
modes = 1;
amps = 1;
ranges = 1;


if ~isempty(varargin) && isscalar(varargin{1})
    amp = varargin{1};
    if length(varargin) > 1 && isscalar(varargin{2})
        range = varargin{2};
    end
end

%% Parse input

for i = 1:length(varargin)
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'com')
        complexKernel = 1;
    end
    
    if strcmpi(varargin{i},'modes') && length(varargin) > i && ~ischar(varargin{i+1})
        modes = varargin{i+1};
    end
    
    if strcmpi(varargin{i},'amps') && length(varargin) > i && ~ischar(varargin{i+1})
        amps = varargin{i+1};
        amps = amps/max(amps); % normalize so that max(amps) = 1
    end
    
    if strcmpi(varargin{i},'ranges') && length(varargin) > i && ~ischar(varargin{i+1})
        ranges = varargin{i+1};
    end

end


if length(modes) > 1
    
    if length(amps) == 1
        amps = repmat(amps,1,length(modes));
    elseif length(amps) ~= length(modes)
        error('modes, amps, and ranges must be same length, or amps and ranges can be scalars that will be copied');
    end
    
    if length(ranges) == 1
        ranges = repmat(ranges,1,length(modes));
    elseif length(ranges) ~= length(modes)
        error('modes, amps, and ranges must be same length, or amps and ranges can be scalars that will be copied');
    end
    
end

%% Generate connectivity kernel(s)

switch n2.fspac
    case 'log'
        F = log2(relFreq(n1.f, n2.f)); % F: log2 of freq ratios
        sigma = abs(log2(range))*ranges;
        per = floor(n1.N/(log2(n1.f(end))-log2(n1.f(1))));
        df = 1/per;
    case 'lin'
        F = diffFreq(n1.f, n2.f); % F: freq differences
        sigma = range*ranges;
        df = n1.df;
end

R = zeros(size(F));
Q = zeros(size(F));

for nn = 1:length(modes)
    switch type
        case 'full'
            R1 = ones(size(F));
            %       R1 = R1/n1.N; % normalize so that sum = 1
        case 'near'
            R1 = amps(nn)*(abs(F-log2(modes(nn))) <= sigma(nn) + eps);
            % JCK had to add eps b/c Matlab is inconsistent at +/- boundaries
            R1 = R1/(2*floor((sigma(nn) + eps)/df)+1); % normalization
        case 'one'
            R1 = eye(n2.N, n1.N);
        case 'gaus'
            R1 = amps(nn)*normpdf(F, log2(modes(nn)), sigma(nn));
            R1 = R1*df; % normalization
            if complexKernel
                Q1 = (pi/2)*(2*normcdf(F, log2(modes(nn)), sigma(nn))-1);
            end
    end
    % R and Q are updated where R1 > R
    ind = find(R1 > R);
    R(ind) = R1(ind);
    if complexKernel && strcmpi(type, 'gaus')
        Q(ind) = Q1(ind);
    end
end

C = R.*exp(1i*Q);
C = amp*C; % Scale to user-specified strength
C(abs(C)<max(max(abs(C)))/100) = 0; % eliminate zeros so we can declare sparse

% =====================================================================================
function RF = relFreq(f1, f2)
% Ratio of f1s (from) in terms of f2s (to)
[F1, F2] = meshgrid(f1, f2);
RF = F1./F2;
% =====================================================================================
function DF = diffFreq(f1, f2)
% f1s (from) minus f2s (to)
[F1, F2] = meshgrid(f1, f2);
DF = F1 - F2;
%% revision history
% 12/15/2011 JCK Combined and modified connectMakeLog and connectMakeLin
%                Added varargin for defining multimodal connections
%                Log, lin fspacs treated equivalently using log2(relFreq), diffFreq