%% example1.m
%
% A one layer network driven with a sinusoidal input. Several parameter
% sets are provided for experimentation with different types of intrinsic
% oscillator dynamics.
%
% https://github.com/MusicDynamicsLab/GrFNNToolbox/wiki/04.-Example-1

%% Choose a parameter set

% alpha =-1; beta1 =  0; beta2 =  0; delta1 = 0; delta2 = 0; neps = 1; % Linear
alpha = 0; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Critical
% alpha = 0; beta1 = -1; beta2 = -1; delta1 = 1; delta2 = 0; neps = 1; % Critical with detuning
% alpha = 1; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Limit Cycle
% alpha =-1; beta1 =  3; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Double limit cycle

%% Make the model
s = stimulusMake(1, 'fcn', [0 50], 40, {'exp'}, [1], .25, 0, ...
    'ramp', 0.02, 1, 'display', 0);

n = networkMake(1, 'hopf', alpha, beta1,  beta2, delta1, delta2, neps, ...
    'log', .125, 8, 601, 'save', 1, 'display', 0, ...
    'Tick', [.5 2/3 3/4 1 4/3 3/2 2], 'zna',.0001);

% n = connectAdd(s, n, 1, 'type', 'allfreq'); % default connection type for stimulus source is '1freq'
n = connectAdd(s, n, 1); % default connection type for stimulus source is '1freq'

M = modelMake(s, n);

options=odeset('RelTol',1e-2);
tic

M.odefun = @ode45;
[T,Z] = M.odefun(@zdotAdaptive,s.ts,n.z,options,M);

% M = M.odefun(M);

toc

M.n{1}.t=T;
M.n{1}.Z=Z.';

%% Display the output
figure; clf; a1 = gca;
figure; clf;
a2 = subplot('Position', [0.08  0.72  0.78 0.22]);
a3 = subplot('Position', [0.08  0.10  0.88 0.50]);

outputDisplay(M, 'net', 1, a1, 'ampx', a2, 'fft', a3, 'oscfft')
