%% example2.m
%
% A simple afferent chain network with no learning
%
% https://github.com/MusicDynamicsLab/GrFNNToolbox/wiki/06.-Example-2

%% Explore different parameter sets
alpha1 = 0.01; beta11 = -1; beta12 =  -10; neps1 = 1; % Layer 1
alpha2 =   -1; beta21 =  4; beta22 =  -3; neps2 = 1; % Layer 2

%% Make the model
s = stimulusMake(1, 'fcn', [0 1], 4000, {'exp'}, [100], .025, 0, 'ramp', 0.01, 1, ...
    'display', 10);

n1 = networkMake(1, 'hopf', alpha1, beta11,  beta12,  0, 0, neps1, ...
    'log', 50, 200, 201, 'save', 1, ...
    'display', 10, 'Tick', [50 67 75 100 133 150 200]);

n2 = networkMake(2, 'hopf', alpha2, beta21,  beta22,  0, 0, neps2, ...
    'log', 50, 200, 201, 'save', 1, ...
    'display', 10, 'Tick', [50 67 75 100 133 150 200]);

n1 = connectAdd(s, n1, 1); % '1freq' connection type by default

C     = connectMake(n1, n2, 'one', 1, 1);
n2    = connectAdd(n1, n2,  C, 'weight', 1, 'type', '1freq');

M = modelMake(s, n1, n2);

tic
M = M.odefun(M);
toc

% %% Display the output
figure(11); clf;
a1 = subplot(2,1,1);
a2 = subplot(2,1,2);

outputDisplay(M,'net',1,a1,'ampx','net',2,a2,'ampx')

figure(12); clf;
a3 = subplot(2,1,1);
a4 = subplot(2,1,2);

outputDisplay(M,'net',1,a3,'fft','net',2,a4,'fft')
