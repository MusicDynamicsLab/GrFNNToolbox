usegpu=1;

%% example1P.m
%
% A one layer network with plastic internal connections (and no input)

%% Network parameters
alpha = 1; beta1 = -1; beta2 = -1000; neps = 1; % Limit Cycle

%% Parameter sets for Hebbian plasiticity
w = .05;
% lambda =  -.1; mu1 =  0; mu2 =  0; ceps =  4; kappa = 1; % Linear learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Critical, stronger nonlinearity
lambda = .001; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Supercritical learning rule

%% Make the model
s = stimulusMake('fcn', [0 100], 40, {'exp'}, [1], [0]);

n = networkMake(1, 'hopf', alpha, beta1,  beta2, 0, 0, neps, ...
    'log', .5, 2, 200, 'channel', 1, 'save', 1, ...
    'display', 0, 'Tick', [.5 .67 .75 1 1.25 1.33 1.50 2]);

n = connectAdd(n, n, [], 'weight', w, 'type', 'all2freq', ...
                         'learn', lambda, mu1, mu2, ceps, kappa, ...
                         'display', 10,'phasedisp', 'save', 500);

if usegpu
    
    M = modelMake(@zdot_gpu, @cdot_gpu, s, n);
    
else
    
    M = modelMake(@zdot, @cdot, s, n);
    
end

%% Run the network

if usegpu
    
    tic
    modelTemp = odeRK4fs_gpu(M);
    toc
    
    M.n{1}.Z=modelTemp.n{1}.Z;
    M.n{1}.con{1}.C3=modelTemp.n{1}.con{1}.C3;
    
else
    
    tic
    M = odeRK4fs(M,s);
    toc
    
end

%% Display the output
if usegpu
    figure(13); clf; a1 = gca;
    figure(14); clf;
else
    figure(11); clf; a1 = gca;
    figure(12); clf;
end
a2 = subplot('Position', [0.08  0.72  0.78 0.22]);
a3 = subplot('Position', [0.08  0.10  0.88 0.50]);

outputDisplay(M, 'net', 1, a1, 'ampx', a2, 'fft', a3, 'oscfft')
