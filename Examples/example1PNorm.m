%% example1PNorm.m
%
% A one layer network with plastic internal connections (and no input)
% A modified version of example1P.m to test buffer normalization

%% Network parameters
% alpha = 1; beta1 = -1; beta2 = -1000; neps = 1; % Limit Cycle
alpha = 1; beta1 = -1; beta2 = -1; neps = 1; % spont amp ~ 0.7

%% Parameter sets for Hebbian plasiticity
w = .05;
% w = .5;
% lambda =  -.1; mu1 =  0; mu2 =  0; ceps =  4; kappa = 1; % Linear learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Critical, stronger nonlinearity
% lambda = .001; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Supercritical learning rule
lambda = .001; mu1 = -1; mu2 = -1; ceps = 1; kappa = .5; % Supercritical learning rule

%% Make the model
s = stimulusMake(1, 'fcn', [0 100], 40, {'exp'}, 1, 0);

n = networkMake(1, 'hopf', alpha, beta1,  beta2, 0, 0, neps, ...
    'log', .5, 2, 201, 'save', 0, ...
    'display', 10, 'Tick', [.5 .67 .75 1 1.25 1.33 1.50 2]);

% unnormalized
n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', ...
    'learn', lambda, mu1, mu2, ceps, kappa, 'display', 10, 'phasedisp', 'save', 500);

% normalized
% n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', 'norm', ...
%     'learn', lambda, mu1, mu2, ceps, kappa, 'display', 10, 'phasedisp', 'save', 500);

% no11, unnormalized
% n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', 'no11', ...
%     'learn', lambda, mu1, mu2, ceps, kappa, 'display', 10, 'phasedisp', 'save', 500);

% no11, normalized
% n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', 'no11', 'norm', ...
%     'learn', lambda, mu1, mu2, ceps, kappa, 'display', 10, 'phasedisp', 'save', 500);

M = modelMake(@zdotNorm, @cdot, s, n);
M.odefun = @odeRK4fsNorm;

tic
M = M.odefun(M);
toc
