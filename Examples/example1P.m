%% example1P.m
%
% A one layer network with plastic internal connections (and no input)
%
% https://github.com/MusicDynamicsLab/GrFNNToolbox/wiki/5.-Example-1-Plastic

%% Network parameters
alpha = 1; beta1 = -1; beta2 = -1000; neps = 1; % Limit Cycle

%% Parameter sets for Hebbian plasiticity
w = .05;
% lambda =  -.1; mu1 =  0; mu2 =  0; ceps =  4; kappa = 1; % Linear learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Critical, stronger nonlinearity
lambda = .001; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Supercritical learning rule

%% Make the model
s = stimulusMake(1, 'fcn', [0 100], 40, {'exp'}, 1, 0);

n = networkMake(1, 'hopf', alpha, beta1,  beta2, 0, 0, neps, ...
    'log', .5, 2, 200, 'save', 1, ...
    'display', 10, 'Tick', [.5 .67 .75 1 1.25 1.33 1.50 2]);

n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', ...
    'learn', lambda, mu1, mu2, ceps, kappa, ...
    'display', 10,'phasedisp', 'save', 500);

M = modelMake(@zdot, @cdot, s, n);
        % The network is not connected to the stimulus, but the model needs
        % a stimulus to get a time vector

%% Run the network
tic
M = odeRK4fs(M);
toc
