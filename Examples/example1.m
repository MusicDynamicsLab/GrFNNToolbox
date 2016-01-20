usegpu=1;

%% example1.m
%
% A one layer network driven with a sinusoidal input. Several parameter
% sets are provided for experimentation with different types of intrinsic
% oscillator dynamics.

%% Choose a parameter set

% alpha =-1; beta1 =  0; beta2 =  0; delta1 = 0; delta2 = 0; neps = 1; % Linear
alpha = 0; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Critical
% alpha = 0; beta1 = -1; beta2 = -1; delta1 = 1; delta2 = 0; neps = 1; % Critical with detuning
% alpha = 1; beta1 = -1; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Limit Cycle
% alpha =-1; beta1 =  3; beta2 = -1; delta1 = 0; delta2 = 0; neps = 1; % Double limit cycle

%% Make the model
s = stimulusMake('fcn', [0 50], 40, {'exp'}, [1], .25, 0, ...
    'ramp', 0.02, 1, 'display', 0);
% stimulusShow(s, 1); drawnow;

n = networkMake(1, 'hopf', alpha, beta1,  beta2, delta1, delta2, neps, ...
                   'log', .5, 2, 200, 'channel', 1, 'save', 1, ...
                   'display', 0, 'Tick', [.5 2/3 3/4 1 4/3 3/2 2]);

               
if usegpu
    
    M = modelMake(@zdot_gpu, @cdot_gpu, s, n);
    
    %% Run the network
    tic
    Mtemp = odeRK4fs_gpu(M);
    toc
    
    M.n{1}.Z=Mtemp.n{1}.Z;
    
else
    
    M = modelMake(@zdot, @cdot, s, n);
    
    %% Run the network
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