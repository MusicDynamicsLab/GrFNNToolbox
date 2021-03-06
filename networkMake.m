%% networkMake
%  n = networkMake(idNumber, varargin)
% 
%  Makes an oscillator network structure for a model
%  idNumber: is a layer id, should match model structure index
%  
%  Attributes can come in any order.
% 
%  Necessary attributes:
%  Model, such as 'hopf'                Next 6 inputs are necessary, parameters of the network: 
%                                       alpha, beta1, beta2, delta1, delta2, epsilon
%  Frequency spacing, 'lin' or 'log'    Next 3 inputs are necessary, min freq of network, 
%                                       max freq of network, and number of oscillators
% 
%  Optional attributes:
%  'scale', 'noScale'                   Use these to specify whether or not to scale parameters
%                                       by natural frequencies. No additional argument.
%  'display'                            Next input is scalar which is the time step interval in 
%                                       between display steps as the network integrates. Default
%                                       is zero.
%  'save'                               Next input is scalar which is the time step interval in 
%                                       between save steps as the network integrates. Default is
%                                       zero.
%  'znaught'                            Next input is a vector of initial conditions for the oscillators 
%                                       overriding the defaults, which are the spontaneous amplitudes
%                                       of the oscillators with a small amount of randomness, and 
%                                       random phase. Can be a scalar if that value should be used 
%                                       for all oscillators.
%  'tick'                               Next input is a vector of axis ticks to be used for this network
%                                       for plotting, if 'display' is nonzero.
% 
%  Example calls:
% 
%   n  = networkMake(3, 'hopf', .1, -10, -1, 0, 0, 1, 'log', 100, 1300, 400, 'display', 20, 'save', 1, 'znaught', z0{3});
%   n  = networkMake(1, 'hopf', 0, -100, -1, 0, 0, .0025, 'log', 20, 20000, 800, 'noScale', 'save', 1, 'channel', 1);
%

%%
function n = networkMake(id, varargin)
%% Set defaults and initialize variables
n.id = id;
n.class = 'network';
n.nClass = 2; % numerical class

n.con = {}; % cell array containing connections targeting this network
n.learnList = []; % indices for learned connections (used in integrator)

models = {'hopf'};              % Can add to this array later
n.model = [];                   % Initialize these to use isempty to error check later
n.fspac = [];

n.nFspac= 0;
n.dStep = 0;                    % Initialize these to zero/empty in case not specified in varargin
n.sStep = 0;
n.tick  = [];
overrideInitialConditions = 0;
userScale = [];

%% Parse input
for i = 1:length(varargin)
    
    if any(strcmpi(varargin{i},models)) && length(varargin) > i + 5 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) ...
            && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5}) && ~ischar(varargin{i+6})
        
        n.model = lower(varargin{i});
        
        switch varargin{i}
            
            case {'hopf'} % 'hopp'? Let's make it 'hopf' ...
                alpha   = varargin{i+1};
                beta1   = varargin{i+2};
                beta2   = varargin{i+3};
                delta1  = varargin{i+4};
                delta2  = varargin{i+5};
                epsilon = varargin{i+6};
                
        end
    end
    
    if ischar(varargin{i}) && any(strcmpi(varargin{i}(1:3),{'lin' 'log'})) && length(varargin) > i + 2 ...
            && isscalar(varargin{i+1}) && isscalar(varargin{i+2}) && isscalar(varargin{i+3})
        
        n.fspac = lower(varargin{i}(1:3));
        lf   = varargin{i+1};            % min freq
        hf   = varargin{i+2};            % max freq
        N    = varargin{i+3};            % number of frequency steps
        n.N  = N;
        
        if strcmpi(n.fspac, 'lin')
            n.nFspac = 1;       % numerical frequency spacing type
            n.scale = 0;
            n.f  = linspace(lf,hf,N)';
            if N > 1
                n.df = abs(n.f(2)-n.f(1)); % to scale for frequency density
            else
                n.df = 1;
            end
        else
            n.nFspac = 2;
            n.scale = 1;
            n.f  = logspace(log10(lf),log10(hf),N)';
            if N > 1
                n.df = abs(log2(n.f(2)/n.f(1)));          % to scale for frequency density
            else
                n.df = 1;
            end
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sca')
        
        userScale = 1;
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'nos')
        
        userScale = 0;
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'dis') && length(varargin) > i && isscalar(varargin{i+1})
        
        n.dStep = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sav') && length(varargin) > i && isscalar(varargin{i+1})
        
        n.sStep = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'zna') && length(varargin) > i && ~ischar(varargin{i+1})
       
        overrideInitialConditions = 1;
        z0 = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'tic') && length(varargin) > i && ~ischar(varargin{i+1})
    
        n.tick = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && ~any(strcmpi(varargin{i},models)) && ~any(strcmpi(varargin{i}(1:3),{'lin' 'log'})) ...
            && ~any(strcmpi(varargin{i}(1:3),{'sca' 'nos'})) && ~strcmpi(varargin{i}(1:3),'dis') ...
            && ~strcmpi(varargin{i}(1:3),'sav') && ~strcmpi(varargin{i}(1:3),'zna') && ~strcmpi(varargin{i}(1:3),'tic')
        
        error(['Unrecognized input to networkMake: ' varargin{i}]);
        
    end
    
end

%% Error check for necessary inputs
if isempty(n.model)
    
    error('Must specify a model and all parameters in networkMake')
    
end


%% Define oscillator parameters
if ~isempty(userScale)
    n.scale = userScale;    % use user-specified scale flag if any
end

switch n.scale
    
    case 0 % no frequency scaling
        n.a  = alpha + 1i*2*pi.*n.f;
        n.b1 = beta1 + 1i*delta1;
        n.b2 = beta2 + 1i*delta2;
        
    case 1 % frequency scaling
        n.a  = (alpha + 1i*2*pi  ).*n.f;  % Redefinition of a, b1 & b2
        n.b1 = (beta1 + 1i*delta1).*n.f;
        n.b2 = (beta2 + 1i*delta2).*n.f;
end

n.e = epsilon;


switch n.model
    case {'hopf', 'hopp'}
        r0 = zeros(size(n.a));
        r = spontAmp(real(n.a(1)), real(n.b1(1)), real(n.b2(1)), n.e);
        r0 = r(end)+r0;
        r0 = r0+.01*randn(size(r0));
        phi0 = 2*pi*randn(size(r0));
        n.z0 = r0.*exp(1i*2*pi*phi0);
        
    otherwise
        error('Unknown model type: %s', n.model);
        
end

if overrideInitialConditions
    if numel(z0) == 1
        n.z0 = repmat(z0,N,1);
    elseif all(size(z0) == [1 N])
        n.z0 = z0.';
    elseif all(size(z0) == [N 1])
        n.z0 = z0;
    else
        error('Length of initial conditions must be 1 or N')
    end        
end

n.z = n.z0;

% Commenting out all former tick stuff and letting matlab decide if not spec'd in varargin
%
% Define ticks and tick labels to be used for plotting
% m   = ceil(n.N/2);                  % middle frequency of network
% switch n.nFspac
%     case 'lin'
%         n.tck = floor(linspace(1, n.N,5));
%     case 'log'
%         per = floor(n.N/(log2(n.f(end))-log2(n.f(1))));
% 
%         tckup = m:per:n.N;
%         tckdn = m:-per:1;
%         n.tck = unique([tckdn,tckup]);
% end
% n.tckl = {};
% for tt = 1:length(n.tck)
%     if n.f(m) < 10
%         n.tckl{tt} = sprintf('%5.1f', n.f(n.tck(tt)));
%     else
%         n.tckl{tt} = sprintf('%d', round(n.f(n.tck(tt))));
%     end
% end

% Older things not to throw away

% model is oscillator-level model-type {'vrd', 'wils', 'hopf'} + model
% parameters
% -for W-C: model = {'wils', a, b, omega, xstar, ystar}
%  parameters c, d of W-C are calcaulated internally based on omega, and
%  the desired fixed point (xstar, ystar)
% -for Hopf: model = {'hopf', a, b1, b2, d1, d2, epsilon} with a, b1, b2,
%  d1, d2, epsilon are all real
% -for VDR = {'vdr', a, b, c, d, A, B}
% freqs is parameters for eigen-frequencies in the oscillator array
% for linear spaced array freqs is a list of {'lin', minimum frequency,
% frequency step-size, and maximum frequency}
% for log-spaced array freqs is a list of {'log', center frequency, number of octaves (on each side),
% oscillators per octave}
%

% switch n.model
%     % Other cases that could be implemented:
%     %   case 'term'
%     %       error('Terman-Wang model not yet implemented');
%     %   case 'fitz'
%     %       error('Fitzhugh-Nagumo model not yet implemented');
%
%     case {'vrd', 'vrd2'}
%         n.a = model{2};
%         n.b = model{3};   % coefficient for the vanderPol term
%         n.c = model{4};   % coefficient for the Rayleigh term
%         n.d = model{5};   % coefficient for the Duffing term
%         n.A = model{6};
%         n.B = model{7};
%     case {'wils', 'wilse'}
%         n.a  = model{2};
%         n.b  = model{3};
%         n.Om = model{4};
%         n.xe = model{5};
%         n.ye = model{6};
%         if length(model)>6
%             n.e = model{7};
%         else
%             n.e = 1.0;
%         end
%
%         [n.c, n.d, n.px, n.py] = wcparams(n.a, n.b, n.Om, n.xe, n.ye);
%         n.mx = n.Om/2/pi;
%         n.my = n.mx;
%
%         n.a1 = -1 + n.a*n.xe*(1-n.xe);
%         n.a2 = -    n.b*n.xe*(1-n.xe);
%         n.a3 =      n.c*n.ye*(1-n.ye);
%         n.a4 = -1 + n.d*n.ye*(1-n.ye);
%
%         n.c1 = 1;
%         n.c3 = 1;
%
% end % end switch/case

% switch n.model
%     case {'hopf', 'hopp'}
%         n.z0 = 1e-10*(1*ones(size(n.f))+i*ones(size(n.f)));
%
%     case {'wils', 'wilse'}
%         x0  = n.xe*ones(size(n.f));
%         y0  = n.ye*ones(size(n.f));
%
%         n.z0 = zeros(2*n.N,1);
%         n.z0(1:2:2*n.N) = x0;
%         n.z0(2:2:2*n.N) = y0;
%     case {'kura'}
%         % start all oscillators at amp=1, phase = 0; (amp=1 is implicit)
%         n.z0 = zeros(size(n.f));
%         n.z0 = zeros(size(n.f)) + n.init_dispersion*(2*pi*rand(size(n.f))-pi);
%
%     otherwise
%         v1  = 1e-10*ones(size(n.f));
%         v2  = 1e-10*ones(size(n.f));
%         n.z0 = zeros(2*n.N,1);
%         n.z0(1:2:2*n.N) = v1;
%         n.z0(2:2:2*n.N) = v2;
% end
