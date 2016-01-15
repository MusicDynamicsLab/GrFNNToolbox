%% modelMake
%  Input dotfunc handle first, then optionally cdot function handle.
%  Then stimulus structure, then each network structure.
%
%  Example calls:
%
%   m = modelMake(@zdot, s, n);
%   m = modelMake(@zdot, @cdot, s, n1, n2, n3);
%
%  collect output network with:
%  n = m.n{1};

%%
function model = modelMake(varargin)

model.dotfunc      = varargin{1};
if isa(varargin{2},'function_handle')
    model.cfun     = varargin{2};
    ind            = 3;
else
    model.cfun     = @cdot;
    ind            = 2;
end

%% Get stimuli first to get time info
stimList = []; % list of stimulus id's
netList = []; % list of network id's
fs = [];
dt = [];
ts = [];
t = [];

for v = ind:length(varargin)
    temp = varargin{v};
    if strcmp(temp.class, 'stim')
        sid = temp.id;
        temp.N = size(temp.x, 1);
        temp.z = temp.x(:,1);   % initialize current state z
        model.s{sid} = temp;
        stimList = [stimList sid];
        
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

model.fs           = fs;
model.dt           = dt;
model.tspan        = ts;
model.t            = t;

%% Get networks and set initial conditions

for v = ind:length(varargin)
    temp = varargin{v};
    if strcmp(temp.class, 'net')
        nid = temp.id;
        model.n{nid} = temp;
        netList = [netList nid];
        
        if temp.sStep > 0
            Nt = ceil(length(t)/temp.sStep);
            model.n{nid}.t = t(1:temp.sStep:length(t));
            model.n{nid}.Z = single(zeros(length(temp.z), Nt));
            model.n{nid}.Z(:,1) = temp.z0;
        else
            model.n{nid}.t = [];
            model.n{nid}.Z = [];
        end
        
        for cx = temp.conLearn
            con = temp.con{cx};
            if con.sStep > 0
                Nt = ceil(length(t)/con.sStep);
                model.n{nid}.con{cx}.t = t(1:con.sStep:length(t));
                model.n{nid}.con{cx}.C3 = single(zeros(size(con.C,1), size(con.C,2), Nt));
                model.n{nid}.con{cx}.C3(:,:,1) = con.C0;
            else
                model.n{nid}.con{cx}.t  = [];
                model.n{nid}.con{cx}.C3 = [];
            end
            
        end
    end
end

model.stimList = sort(stimList);
model.netList = sort(netList);

%% Check if all connections are valid and at least one network gets stimulus
stimcount = 0;
for j = model.netList
    net = model.n{j};
    if net.ext   % connect first stimulus if ext is nonzero (backward compatible)
        stim1 = model.s{model.stimList(1)};
        if stim1.N == 1 % single channel input
            C = ones(net.N, 1);
        else    % multichannel input
            C = zeros(net.N, stim1.N);
            C(sub2ind(size(C), 1:net.N, net.ext)) = 1;
        end
        model.n{j} = connectAdd(stim1, net, C, 'weight', 1, 'type', stim1.inputType);
    end
    
    for k = 1:length(model.n{j}.con)
        con = model.n{j}.con{k};
        if strcmp(con.sourceClass, 'stim')
            if ~ismember(con.source, model.stimList)
                error(['Input to Network ' num2str(j) ' is missing (Stimulus ' num2str(k) ')'])
            end
            stimcount = stimcount + 1;
        end
        if strcmp(con.sourceClass, 'net')
            if ~ismember(con.source, model.netList)
                error(['Input to Network ' num2str(j) ' is missing (Network ' num2str(k) ')'])
            end
        end
    end
end

if ~stimcount   % if no network is connected to stimulus
    stim1 = model.s{model.stimList(1)};
    net1 = model.n{model.netList(1)};
    if stim1.N == 1 % single channel input
        C = ones(net1.N, 1);
    else    % multichannel input
        C = eye(net1.N, stim1.N);
    end
    model.n{model.netList(1)} = connectAdd(stim1, net1, C, 'weight', 1, 'type', '1freq');
    disp({'At least one network must be connected to a stimulus.',['modelMake connected Network ' num2str(model.netList(1)) ' to Stimulus ' num2str(model.stimList(1)) '.']})
end

%% Cast everything as complex and single

for nx = model.netList
    model.n{nx}.z0 = castCS(model.n{nx}.z0);
    model.n{nx}.z  = castCS(model.n{nx}.z);
    model.n{nx}.Z  = castCS(model.n{nx}.Z);
    model.n{nx}.a  = castCS(model.n{nx}.a);
    model.n{nx}.b1 = castCS(model.n{nx}.b1);
    model.n{nx}.b2 = castCS(model.n{nx}.b2);
    model.n{nx}.e  = castCS(model.n{nx}.e);
    for cx = 1:length(model.n{nx}.con)
        if any(model.n{nx}.conLearn) && any(model.n{nx}.conLearn == cx)
            model.n{nx}.con{cx}.C0 = castCS(model.n{nx}.con{cx}.C0);
            model.n{nx}.con{cx}.C  = castCS(model.n{nx}.con{cx}.C);
            model.n{nx}.con{cx}.C3 = castCS(model.n{nx}.con{cx}.C3);
        end
        model.n{nx}.con{cx}.w      = castCS(model.n{nx}.con{cx}.w);
        model.n{nx}.con{cx}.lambda = castCS(model.n{nx}.con{cx}.lambda);
        model.n{nx}.con{cx}.mu1    = castCS(model.n{nx}.con{cx}.mu1);
        model.n{nx}.con{cx}.mu2    = castCS(model.n{nx}.con{cx}.mu2);
        model.n{nx}.con{cx}.kappa  = castCS(model.n{nx}.con{cx}.kappa);
        model.n{nx}.con{cx}.e      = castCS(model.n{nx}.con{cx}.e);
    end
end


%% If dotfunc (override) option is empty, then use base dotfunc from network and oscillator-model

if isempty(model.dotfunc)
    n = varargin{1}; % for now, restrict to only one osc-model
    if strcmp(n.model, 'vdp')
        model.dotfunc = @zdotv;
    end
    if strcmp(n.model, 'wc')
        model.dotfunc = @zdotw_sc;
    end
    if strcmp(n.model, 'wce')
        model.dotfunc = @zdotw_sc;
    end
    if strcmp(n.model, 'hopft')
        model.dotfunc = @zdotw_sc;
    end
    if strcmp(n.model, 'hopfx')
        model.dotfunc = @zdotw_sc;
    end
    
end
