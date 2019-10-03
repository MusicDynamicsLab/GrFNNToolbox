%% example1.m
%
% A variant of example1.m that tests buffer normalization.

%% Choose a parameter set
% alpha =-1; beta1 =  0; beta2 =  0; delta1 = 0; delta2 = 0; neps = 1; % Linear
alpha = 0; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Critical
% alpha = 0; beta1 = -1; beta2 = -1; delta1 = 1; delta2 = 0; neps = 1; % Critical with detuning
% alpha = 1; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Limit Cycle
% alpha =-1; beta1 =  3; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Double limit cycle

%% Make the model
s = stimulusMake(1, 'fcn', [0 50], 40, {'exp'}, [1], .7, 0, ...
    'ramp', 0.02, 1, 'display', 0);

n = networkMake(1, 'hopf', alpha, beta1,  beta2, delta1, delta2, neps, ...
    'log', .5, 2, 201, 'save', 1, 'display', 10, ...
    'Tick', [.5 2/3 3/4 1 4/3 3/2 2]);

n = connectAdd(s, n, 1, 'type', 'all2freq'); % unnormalized
% n = connectAdd(s, n, 1, 'type', 'all2freq', 'norm'); % normalized
% n = connectAdd(s, n, 1, 'type', 'all2freq', 'no11'); % no 11, unnormalized
% n = connectAdd(s, n, 1, 'type', 'all2freq', 'norm', 'no11'); % no11, normalize

M = modelMake(@zdotNorm, s, n);
M.odefun = @odeRK4fsNorm;

tic
M = M.odefun(M);
toc

%% Display the output
figure(11); clf; a1 = gca;
figure(12); clf;
a2 = subplot('Position', [0.08  0.72  0.78 0.22]);
a3 = subplot('Position', [0.08  0.10  0.88 0.50]);

outputDisplay(M, 'net', 1, a1, 'ampx', a2, 'fft', a3, 'oscfft')
