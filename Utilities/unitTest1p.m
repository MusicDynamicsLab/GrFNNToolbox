%% unitTest1p.m
%
% A one layer network driven with a sinusoidal input. Several parameter
% sets are provided for experimentation with different types of intrinsic
% oscillator dynamics.
% 
% This script is similar to example1p.m from the GrFNN Toolbox, but instead
% compares the results from the CPU and GPU implementations.

connectionTypes = {'1freq', 'all2freq', 'allfreq', 'active'};
% Need to have a separate array for use in the output table since
% Matlab uses these as "variable" names which cannot start with "1".
connectionNames = {'onefreq', 'all2freq', 'allfreq', 'active'};

aLin  = -1; aCrit  = 0; aCritDetune =  0; aLC = 1; aDLC = -1;
b1Lin =  0; b1Crit =-1; b1CritDetune= -1; b1LC=-1; b1DLC = 3;
b2Lin =  0; b2Crit =-1; b2CritDetune= -1; b2LC=-1e3; b2DLC =-1;
d1Lin =  0; d1Crit = 0; d1CritDetune=  1; d1LC= 0; d1DLC = 0;
d2 = 0;
eLin  =  1; eCrit  = 1; eCritDetune =  1; eLC = 1; eDLC  = 1;

params = struct('alpha', [aLin, aCrit, aCritDetune, aLC, aDLC],...
    'beta1', [b1Lin, b1Crit, b1CritDetune, b1LC, b1DLC],...
    'beta2', [b2Lin, b2Crit, b2CritDetune, b2LC, b2DLC],...
    'delta1', [d1Lin, d1Crit, d1CritDetune, d1LC, d1DLC],...
    'delta2', [d2, d2, d2, d2, d2],...
    'eps', [eLin, eCrit, eCritDetune, eLC, eDLC]);
pr = 4; % choose which parameter regime to use

% Parameters for Hebbian plasticity
lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule

w = 0; % Connection weight

fs = 40;

numConnectionTypes = length(connectionTypes);

nansPresent1to1      = zeros(1, numConnectionTypes);
nansPresentno11      = zeros(1, numConnectionTypes);

runTimes        = zeros(2, numConnectionTypes * 2);

for ct = 1:numConnectionTypes      % loop over connection types
    disp(['Connection type: ' connectionTypes{ct}]);
    for oto = 0:1 % run with and without one to one connections
        
        %% Make the model
        s = stimulusMake(1, 'fcn', [0 50], fs, {'exp'}, [1], 0, 0, ...
            'ramp', 0.02, 1, 'display', 10);
        
        n = networkMake(1, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
            params.delta1(pr), params.delta2(pr), params.eps(pr),...
            'log', .5, 2, 200, 'save', 1, ...
            'display', 10, 'Tick', [.5 2/3 3/4 1 4/3 3/2 2]);
        
%         C = ones(n.N, s.N);
%         n = connectAdd(s, n, C, 'type', '1freq');
        
        if ~oto
            n = connectAdd(n, n, [], 'type', connectionTypes{ct}, 'no11',...
                'learn', lambda, mu1, mu2, ceps, kappa, 'weight', w, 'display', 100);
        else
            n = connectAdd(n, n, [], 'type', connectionTypes{ct},...
                'learn', lambda, mu1, mu2, ceps, kappa, 'weight', w, 'display', 100);
        end
        
        %
        % Run the network (cpu version)
        %
    
        Mcpu = modelMake(@zdot, @cdot, s, n);
        
        tic
        Mcpu = odeRK4fs(Mcpu);
        runTimes(1, ct + oto) = toc;
      
        Z = Mcpu.n{1}.Z;
        
        
        %
        % Process results
        %
        
        
        if oto
            nansPresent1to1(ct)      = any(isnan(Z(:)));      
        else
            nansPresentno11(ct)      = any(isnan(Z(:)));
        end
    end
end

%
% Display results
%

cpuRuntimeTotal = sum(runTimes(1, :));
cpuRuntimeAvg   = mean(runTimes(1, :));

disp(' ');
disp(sprintf('CPU runtime %0.2fs mean %0.2fs', cpuRuntimeTotal, cpuRuntimeAvg));

output = [nansPresent1to1; nansPresentno11];
metricLabels = {'NaNs Present (oto)','NaNs Present (no11)'};

T = array2table(output,...
    'VariableNames',connectionNames,...
    'RowNames',metricLabels);

disp(' ');
disp(T);
