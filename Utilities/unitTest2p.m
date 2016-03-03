%% unitTest2p.m
%
% Unit test script for two layer network configurations with learning.
%

%% Parameters
connectionTypes = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};

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
pr = 2; % choose which parameter regime to use

% Parameters for Hebbian plasticity
lambdaLin =  -.1; mu1Lin =  0; mu2Lin =  0; cepsLin =  4; kappaLin = 1; % Linear learning rule
lambdaC =   0; mu1C = -1; mu2C = -50; cepsC =  4; kappaC = 1; % Critical learning rule
lambdaSC = .001; mu1SC = -1; mu2SC = -50; cepsSC = 16; kappaSC = 1; % Supercritical learning rule
learningParams = struct('lambda', [lambdaLin, lambdaC, lambdaSC],...
    'mu1', [mu1Lin, mu1C, mu1SC],...
    'mu2', [mu2Lin, mu2C, mu2SC],...
    'ceps', [cepsLin, cepsC, cepsSC],...
    'kappa', [kappaLin, kappaC, kappaSC]);

Fs = 160;
count = 0;
dispRate = 10;
w = 0.05; % Connection weight

for ct = 1:length(connectionTypes)
    disp(['---Connection type: ' connectionTypes{ct} '---']);
    for oto = 0:1 % run with and without one to one connections
        switch oto
            case 0
                disp('no 1to1 connections')
            case 1
                disp('including 1to1 connections')
        end
        for lr = 1:3 % loop over learning regimes
            switch lr
                case 1
                    fprintf('- Linear learning regime running...')
                case 2
                    fprintf('- Critical learning regime running...')
                case 3
                    fprintf('- Supercritical learning regime running...')
            end
            
            count = count+1;
            
            %% Make the model
            s = stimulusMake(1, 'fcn', [0 10], Fs, {'exp'}, [2], .025, 0,...
                'ramp', 0.01, 1, 'display', dispRate);
            
            n1 = networkMake(1, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 4, 200, 'save', 1, ...
                'display', dispRate);
            n2 = networkMake(2, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 8, 100, 'save', 1, ...
                'display', dispRate);
            
            C = connectMake(n1, n2, 'one', 1, 1);
            
            if oto
                n1 = connectAdd(s, n1, 1);
                n2 = connectAdd(n2, n2, [], 'weight', w, 'type', '2freq', ...
                    'learn', learningParams.lambda(lr), learningParams.mu1(lr),...
                    learningParams.mu2(lr), learningParams.ceps(lr),...
                    learningParams.kappa(lr), ...
                    'display', dispRate, 'phasedisp', 'save', 1000);
            else
                n1 = connectAdd(s, n1, 1, 'no11');
                n2 = connectAdd(n2, n2, [], 'weight', w, 'type', '2freq', 'no11', ...
                    'learn', learningParams.lambda(lr), learningParams.mu1(lr),...
                    learningParams.mu2(lr), learningParams.ceps(lr),...
                    learningParams.kappa(lr), ...
                    'display', dispRate, 'phasedisp', 'save', 1000);
            end
            
            M = modelMake(@zdot, @cdot, s, n1, n2);
            
            %% Run the network
            tic
            M = odeRK4fs(M);
            runTimes(count) = toc;
            Z = M.n{1}.Z;
            
            if oto
                nansPresent1to1(ct) = any(isnan(Z(:)));
                if (any(isnan(Z(:))))
                    disp('Warning, NaNs present');
                else
                    disp('OK, no NaNs');
                end
            else
                nansPresentno11(ct) = any(isnan(Z(:)));
                if (any(isnan(Z(:))))
                    disp('Warning, NaNs present');
                else
                    disp('OK, no NaNs');
                end
            end
        end
    end
end

%
% Display results
%

if sum(nansPresentno11)>0 || sum(nansPresent1to1)>0
    disp(' ')
    disp('NaNs present - check above output results')
else
    disp(' ')
    disp('No NaNs present in any of the above simulations')
end