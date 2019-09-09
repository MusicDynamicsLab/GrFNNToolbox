%% modelMake
%  M = modelMake(varargin)
%
%  Creates a model, M, out of stimulus and networks in varargin. Network
%  indices in the created model follow network id's assigned in networkMake
%  calls. Function handles for network and connection dot functions can be
%  specified as first two input arguments (zdot.m and cdot.m are default,
%  and you have to specify network dot function to specify connection
%  dot function).
%
%  Example calls:
%
%   M = modelMake(s, n1, n2);
%   M = modelMake(@zdot, s, n);
%   M = modelMake(@zdot, @cdot, s, n1, n2, n3);
%

%%
function model = modelMake(varargin)

if isa(varargin{1},'function_handle')
    model.zfun = varargin{1};
    if isa(varargin{2},'function_handle')
        model.cfun = varargin{2};
        ind = 3;
    else
        model.cfun = @cdot;
        ind = 2;
    end
else
    model.zfun = @zdot;
    model.cfun = @cdot;
    ind = 1;
end

model.odefun = @odeRK4fs;   % default ode function
model.usingGpu = 0;

%% Get stimuli first to get time info
stimListAll = []; % list of all stimulus id's
stimList = []; % list of stimulus id's, only those used as source
netList = []; % list of network id's
fs = [];
dt = [];
ts = [];
t = [];

for v = ind:length(varargin)
    temp = varargin{v};
    
    if ischar(temp) && strcmpi(temp, 'usegpu')
        model.odefun = @odeRK4fs_gpu;
		model.usingGpu = 1;
    end
    
    if ~ischar(temp) && temp.nClass == 1
        sid = temp.id;
        temp.N = size(temp.x, 1);   % in case extra channel was added
        temp.z = temp.x(:,1);   % initialize current state z
        model.s{sid} = temp;
        stimListAll = [stimListAll sid];
        
        if isempty(fs)
            fs = temp.fs;
            dt = temp.dt;
            ts = temp.ts;
            t = temp.t;
        else
            if ~isequal(fs, temp.fs) || ~isequal(ts, temp.ts)
                error('Time vectors for stimuli must be identical')
            end
        end
    end
end

if isempty(stimListAll)   % if no stimulus is given
    error('No stimulus is given.')
end

model.fs           = fs;
model.dt           = dt;
model.t            = t;
model.iSpan        = [1 length(t)];
model.iStep        = 1;

%% Get networks and set initial conditions

for v = ind:length(varargin)
    temp = varargin{v};
    if ~ischar(temp) && temp.nClass == 2
        nid = temp.id;
        model.n{nid} = temp;
        netList = [netList nid];
        
        if temp.sStep > 0
            Nt = ceil(length(t)/temp.sStep);
            model.n{nid}.t = t(1:temp.sStep:length(t));
            model.n{nid}.Z = zeros(length(temp.z), Nt);
            model.n{nid}.Z(:,1) = temp.z0;
        else
            model.n{nid}.t = [];
            model.n{nid}.Z = [];
        end
        
        for cx = temp.learnList
            con = temp.con{cx};
            if con.sStep > 0
                Nt = ceil(length(t)/con.sStep);
                model.n{nid}.con{cx}.t = t(1:con.sStep:length(t));
                model.n{nid}.con{cx}.C3 = zeros(size(con.C,1), size(con.C,2), Nt);
                model.n{nid}.con{cx}.C3(:,:,1) = con.C0;
            else
                model.n{nid}.con{cx}.t  = [];
                model.n{nid}.con{cx}.C3 = [];
            end
            
        end
    end
end

netList = sort(netList);
stimListAll = sort(stimListAll);

%% Check if all connections are valid and at least one network gets stimulus
for j = netList
    n = model.n{j};
    for k = 1:length(n.con)
        con = n.con{k};
        if con.nSourceClass == 1
            if ~ismember(con.source, stimListAll)
                error(['Input to Network ' num2str(j) ' is missing (Stimulus ' num2str(k) ')'])
            else
                stimList = [stimList con.source];
            end
        end
        if con.nSourceClass == 2
            if ~ismember(con.source, netList)
                error(['Input to Network ' num2str(j) ' is missing (Network ' num2str(k) ')'])
            end
        end
        % Initialize normalization buffer
        if con.norm
            Lbuf = round(fs/min(n.f))*4; % buffer length
            con.Lbuf = Lbuf; % buffer size
            con.buffer = single(zeros(1, Lbuf));
            con.hx = 1; % head position
            beta2 = n.b2./(n.f.^n.scale);
            beta2 = max(beta2);
            con.thrNorm = -beta2/(2*n.e^1.5); % threshold for normalization
            model.n{j}.con{k} = con;
        end
    end
end

model.stimListAll = stimListAll;
model.stimList = sort(stimList);
model.netList = netList;

if isempty(stimList)   % if no network is connected to stimulus
    disp('Warning (modelMake): No network is connected to stimulus.')
end

%% Cast everything as complex and single
model.iSpan    = single(model.iSpan);
model.iStep    = single(model.iStep);
model.dt       = single(model.dt);
for sx = model.stimList
    model.s{sx}.x = castCS(model.s{sx}.x);
end
for nx = model.netList
    model.n{nx}.z0 = castCS(model.n{nx}.z0);
    model.n{nx}.z  = castCS(model.n{nx}.z);
    model.n{nx}.Z  = castCS(model.n{nx}.Z);
    model.n{nx}.a  = castCS(model.n{nx}.a);
    model.n{nx}.b1 = castCS(model.n{nx}.b1);
    model.n{nx}.b2 = castCS(model.n{nx}.b2);
    model.n{nx}.e  = castCS(model.n{nx}.e);
    for cx = 1:length(model.n{nx}.con)
        if model.n{nx}.con{cx}.learn
            model.n{nx}.con{cx}.C3 = castCS(model.n{nx}.con{cx}.C3);
		end
		% DMR Do not castCS on NUM, NUM1, NUM2, and DEN. This will crash the C code (mex).
		% But we _do_ want singles.
        if isfield(model.n{nx}.con{cx},'NUM')
            model.n{nx}.con{cx}.NUM = single(model.n{nx}.con{cx}.NUM);
            model.n{nx}.con{cx}.DEN = single(model.n{nx}.con{cx}.DEN);
        end
        if isfield(model.n{nx}.con{cx},'NUM1')
            model.n{nx}.con{cx}.NUM1 = single(model.n{nx}.con{cx}.NUM1);
            model.n{nx}.con{cx}.NUM2 = single(model.n{nx}.con{cx}.NUM2);
            model.n{nx}.con{cx}.DEN  = single(model.n{nx}.con{cx}.DEN);
        end
        model.n{nx}.con{cx}.C      = castCS(model.n{nx}.con{cx}.C);
        model.n{nx}.con{cx}.C0     = castCS(model.n{nx}.con{cx}.C0);
        model.n{nx}.con{cx}.w      = castCS(model.n{nx}.con{cx}.w);
        model.n{nx}.con{cx}.lambda = castCS(model.n{nx}.con{cx}.lambda);
        model.n{nx}.con{cx}.mu1    = castCS(model.n{nx}.con{cx}.mu1);
        model.n{nx}.con{cx}.mu2    = castCS(model.n{nx}.con{cx}.mu2);
        model.n{nx}.con{cx}.kappa  = castCS(model.n{nx}.con{cx}.kappa);
        model.n{nx}.con{cx}.e      = castCS(model.n{nx}.con{cx}.e);
    end
end


%% If dotfunc (override) option is empty, then use base dotfunc from network and oscillator-model

if isempty(model.zfun)
    n = varargin{1}; % for now, restrict to only one osc-model
    if strcmp(n.model, 'vdp')
        model.zfun = @zdotv;
    end
    if strcmp(n.model, 'wc')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'wce')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'hopft')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'hopfx')
        model.zfun = @zdotw_sc;
    end
    
end
