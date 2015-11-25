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
    s              = varargin{3};
    ind            = 4;
else
    model.cfun     = @cdot;
    s              = varargin{2};
    ind            = 3;
end
model.fs           = s.fs;
model.dt           = s.dt;
model.tspan        = s.ts;

%% Make initial conditions. The varargs are the networks.
model.netList = []; % list of network id's

for v = ind:length(varargin)
    
    nid = varargin{v}.id;
    model.n{nid} = varargin{v};
    model.netList = [model.netList nid];
    
    t = s.t;
    
    if ~isempty(t) && model.n{nid}.sStep > 0
        Nt = ceil(length(t)/model.n{nid}.sStep);
        model.n{nid}.t = t(1:model.n{nid}.sStep:length(t));
        model.n{nid}.Z = single(zeros(length(model.n{nid}.z), Nt));
        model.n{nid}.Z(:,1) = model.n{nid}.z0;
    else
        model.n{nid}.t = [];
        model.n{nid}.Z = [];
    end
    
    for cx = model.n{nid}.conLearn
        
        if ~isempty(t) && model.n{nid}.con{cx}.sStep > 0
            Nt = ceil(length(t)/model.n{nid}.con{cx}.sStep);
            model.n{nid}.con{cx}.t = t(1:model.n{nid}.con{cx}.sStep:length(t));
            model.n{nid}.con{cx}.C3 = single(zeros(size(model.n{nid}.con{cx}.C,1), size(model.n{nid}.con{cx}.C,2), Nt));
            model.n{nid}.con{cx}.C3(:,:,1) = model.n{nid}.con{cx}.C0;
        else
            model.n{nid}.con{cx}.t  = [];
            model.n{nid}.con{cx}.C3 = [];
        end
        
    end
    
end
model.netList = sort(model.netList);

% Roll thru networks and make sure at least one is connected to stimulus
stimcount = 0;
for j = model.netList
    stimcount = stimcount + model.n{j}.ext;
end
if ~stimcount
    model.n{model.netList(1)}.ext = 1;
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
