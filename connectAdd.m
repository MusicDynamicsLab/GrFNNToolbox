%% connectAdd
%  function n = connectAdd(n1, n2, C, varargin)
% 
%  Updates and returns netTo, with new connectivity from netFrom, with
%  connectivity matrix specified by pattern.
%
%  n2 is destination, i.e. 'to'
%
%  Connection struct:
%    conn.source    : where info is coming from
%    conn.type      : afferent, efferent, or internal
%    conn.learn     : whether or not to learn this matrix
%    conn.C   : connectivity matrix
%    conn.lambda, mu1, mu2, kappa, e1 e2: learning rule parameters, can be scalar or matrix
%                          except e1 & e2, which must be scalar
% 
%  First three inputs are required, all others are optional and can come in any order.
%  Attribute 'type' takes one arguemnt after it, which is connection type.
%  Options for 'type' are '1freq', '2freq', '3freq', '3freqAll', 'All2freq', and
%  'Allfreq'. 'All2freq' is the full nonlinearity with single-factor P(),
%  'Allfreq' is the full nonlinearity with two-factor P().
%  Attribute 'learn' takes five arguments after it, the parameters for the learning rule: lambda, mu1, mu2, epsilon, and kappa
%  Attribute 'weight' takes one argument after it, a weight multiplied to the connectivity matrix
%  Attribute 'display' takes one argument after it, the time step interval at which to display the connectivity matrix during integration. Default is zero
%  Attribute 'save' takes one argument after it, the time step interval at which to save the connectivity matrix in the global model variable. Default is zero
%  Single input 'no11' subtracts out all 1-to-1 and subsequent n-to-n
%  monomials from the cdot equation.
%
%  Example calls:
%
%   n2 = connectAdd(n1, n2, con, 'learn', 0, -100000000, -1, .5, .1, 'display', 20, 'type', '1freq');
%   n2 = connectAdd(n1, n2, bm2oc, 'learn', 0, 0, 0, 0, 0);
% 
%

%%
function n = connectAdd(n1, n2, C, varargin)

con.source = n1.id;

types = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};         % Can add to this array later; default to 'Allfreq' for now
con.type = 'Allfreq';

w         = 1;           % Initialize these in case not specified in varargin
con.dStep = 0;
con.sStep = 0;

%% Initialize learning parameters to zero

lambda   = 0; 
mu1      = 0; 
mu2      = 0; 
epsilon  = 0; 
kappa    = 0;
con.no11 = 0;
con.tol  = .01;
con.phaseDisp = 0;

%% Parse input

for i = 1:length(varargin)
    
    if strcmpi(varargin{i},'learn') && length(varargin) > i + 4 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5})
        
        lambda   = varargin{i+1};
        mu1      = varargin{i+2};
        mu2      = varargin{i+3};
        epsilon  = varargin{i+4};
        kappa    = varargin{i+5};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'typ') && length(varargin) > i && ischar(varargin{i+1}) && any(strcmpi(varargin{i+1},types))
        
        con.type = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'2fr') && length(varargin) > i && isscalar(varargin{i+1})
        
        con.tol = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'wei') && length(varargin) > i && isnumeric(varargin{i+1})
        
        w = varargin{i+1};
        
    end
    
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'dis') && length(varargin) > i && isscalar(varargin{i+1})
        
        con.dStep = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sav') && length(varargin) > i && isscalar(varargin{i+1})
        
        con.sStep = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'no11')
        
        con.no11 = 1;
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'phasedisp')
        
        con.phaseDisp = 1;
        
    end
    
    if ischar(varargin{i}) && ~strcmpi(varargin{i},'learn') && ~strcmpi(varargin{i}(1:3),'typ') && ~strcmpi(varargin{i}(1:3),'wei') && ~strcmpi(varargin{i}(1:3),'dis') && ~strcmpi(varargin{i}(1:3),'sav') && ~any(strcmpi(varargin{i},types)) && ~strcmpi(varargin{i},'no11') && ~strcmpi(varargin{i},'phasedisp')
        
        error(['Unrecognized input to connectAdd: ' varargin{i}]) 
        
    end
    
end

%% Connection types

F1 = repmat(n1.f', n2.N, 1);
F2 = repmat(n2.f, 1, n1.N);

switch lower(con.type)
    
    case '1freq' % one-frequency monomials
        F = (F1 + F2)/2;
        
    case '2freq' % two-frequency monomials
        R = F2./F1;
        sz = size(R);
        [N, D] = fareyratio(R(:)', con.tol);
        N = reshape(N,sz);
        D = reshape(D,sz);
        con.NUM = N;
        con.DEN = D;
        F = (N.*F1 + D.*F2)./(N + D);
        
    case '3freq' % three-frequency monomials ALL MONONMIALS, UP TO SPECIFIED ORDER
        if n1.id == n2.id
           [X1i, X2i, Zi, N1, N2, D, CON1, CON2, mask] = threeFreqMatsAll(n1.f);
        else
           [X1i, X2i, Zi, N1, N2, D, CON1, CON2, mask] = threeFreqMatsAll(n1.f, n2.f);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = CON1;
        con.CON2 = CON2;
        con.NUM1 = N1;
        con.NUM2 = N2;
        con.DEN = D;
        con.mask = mask;
        F = (abs(N1).*n1.f(X1i) + abs(N2).*n1.f(X2i) + D.*n2.f(Zi)) ...
            ./(abs(N1) + abs(N2) + D);
        
    case '3freqall' % three-frequency monomials ALL FREQUENCIES, LOWEST ORDER MONOMIAL
        if n1.id == n2.id
            [X1i, X2i, Zi] = inputShapeInternal(n1.N);
            [N1,  N2,  D] = inputExponentsInternal(n1.f, X1i, X2i, Zi);
        else
            [X1i, X2i, Zi] = inputShapeOther(n1.N, n2.N);
            [N1,  N2,  D] = inputExponentsOther(n1.f, n2.f, X1i, X2i, Zi);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = (N1<0);
        con.CON2 = (N2<0);
        con.NUM1 = abs(N1);
        con.NUM2 = abs(N2);
        con.DEN = D;
        F = (abs(N1).*n1.f(X1i) + abs(N2).*n1.f(X2i) + D.*n2.f(Zi)) ...
            ./(abs(N1) + abs(N2) + D);
        
    case {'all2freq','allfreq','active'} % Full series of resonant monomials
        F = (2*F1.*F2)./(F1+F2);
end

%% Scaling for log networks

if strcmpi(n2.fspac,'log')
    con.F = F;
else
    con.F = ones(size(F));
end

%% Set learning parameters

if any(kappa) % JCK: replaces if sum(sum(kappa))
    con.learn = 1;
else
    con.learn = 0;
end

if strcmpi(n2.fspac,'log')
    con.w = w.*n2.f;
    con.lambda = lambda.*F;
    con.mu1 = mu1.*F;
    con.mu2 = mu2.*F;
    con.kappa = kappa.*F;
else
    con.w = w;
    con.lambda = lambda;
    con.mu1 = mu1;
    con.mu2 = mu2;
    con.kappa = kappa;
end
con.e = epsilon;

%% Connection State and Initial Conditions
%       C0: initial conditions
%        C: connection matrix (instantaneous)
%        t: times saved in 3D memory matrix
%       C3: state memory (3D matrix: frequency x frequency x time)

if isempty(C)
    if strcmpi(con.type(1:5), '3freq')
        A0 = zeros(size(con.NUM1));
    else
        A0 = zeros(size(con.F));
    end
    A = spontAmp(real(con.lambda(1,1)), real(con.mu1(1,1)), real(con.mu2(1,1)), con.e);
    A0 = A0+min(A);
    A0 = A0.*(1 +.01*randn(size(A0)));
    theta0 = randn(size(A0));
    con.C0 = A0.*exp(1i*2*pi*theta0);
else
    con.C0 = C;
end

con.C = con.C0;


%% Input
%      n1: index into n1 (the network in the signal is coming from)
%      n2: index into n2 (the network in the signal is going to)
con.n1 = n1.id;
con.n2 = n2.id;


%% Masking
%      This process operates on the initial condition and parameter
%      matrices to hange the way things work (i.e.
%
%      mask = ~eye(size(C));          % Don't learn self-connections
%                                       This could be inverted guassian
%                                       too

mask = 1;
if strcmpi(con.type, '2freq') && con.no11
    mask = ~(con.NUM==1 & con.DEN==1);
elseif strcmpi(con.type, '3freq')
    mask = con.mask;
elseif con.n1 == con.n2
    mask = ~eye(size(con.C));          % ... Won't learn self-connections
end

%% Return network n2
n = n2;
con.lambda = con.lambda .* mask;
con.mu1 = con.mu1 .* mask;
con.mu2 = con.mu2 .* mask;
con.kappa = con.kappa .* mask;
con.C0    = con.C0    .* mask;
con.C     = con.C     .* mask;

n.con{length(n.con)+1} = con;

%% Add index to conLearn if learn
if con.learn
    n.conLearn = [n.conLearn length(n.con)];
end
