%% unitTest2p.m
%
% Unit test script for two layer network configurations with learning.
%

connectionTypes = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};

aLin  = -1; aCrit  = 0; aCritDetune =  0; aLC = 1;   aDLC = -1;
b1Lin =  0; b1Crit =-1; b1CritDetune= -1; b1LC=-1;   b1DLC = 3;
b2Lin =  0; b2Crit =-1; b2CritDetune= -1; b2LC=-1e3; b2DLC =-1;
d1Lin =  0; d1Crit = 0; d1CritDetune=  1; d1LC= 0;   d1DLC = 0;
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
lambdaLin =  -.1; mu1Lin =  0; mu2Lin =  0; cepsLin =  4; kappaLin = 1; % Linear learning rule
lambdaC =   0; mu1C = -1; mu2C = -50; cepsC =  4; kappaC = 1; % Critical learning rule
lambdaSC = .001; mu1SC = -1; mu2SC = -50; cepsSC = 16; kappaSC = 1; % Supercritical learning rule
learningParams = struct('lambda', [lambdaLin, lambdaC, lambdaSC],...
    'mu1', [mu1Lin, mu1C, mu1SC],...
    'mu2', [mu2Lin, mu2C, mu2SC],...
    'ceps', [cepsLin, cepsC, cepsSC],...
    'kappa', [kappaLin, kappaC, kappaSC]);

Fs = 160;
dispRate = 10;
w = 0.05; % Connection weight

for ct = 5%1:length(connectionTypes)
    disp(' ');
    disp(['---Connection type: ' connectionTypes{ct} '---']);
    for oto = 1:-1:0 % run with and without one to one connections
        if oto==1
            disp(' ')
            disp('including 1:1 connections')
        elseif ~oto && (ct==1 || ct==3)
            disp(' ')
            disp('testing no11 not necessary for this connection type')
            continue;
        else
            disp(' ')
            disp('no 1:1 connections')
        end
        for lr = 3%1:3 % loop over learning regimes
            switch lr
                case 1
                    fprintf('- Linear learning regime running...')
                case 2
                    fprintf('- Critical learning regime running...')
                case 3
                    fprintf('- Supercritical learning regime running...')
            end
            
            Nosc = 201;
            if ct==4 % if 3freqall
                Nosc = (Nosc-1)/4 + 1;
            end
            
            %% Make the model
            s = stimulusMake(1, 'fcn', [0 10], Fs, {'exp'}, [2], .25, 0,...
                'ramp', 0.01, 1, 'display', dispRate);
            
            n1 = networkMake(1, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 8, Nosc, 'save', 1, 'Tick', [.5 1 2 4 8],...
                'display', dispRate);
            n2 = networkMake(2, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 8, (Nosc-1)/2+1, 'save', 1, 'Tick', [.5 1 2 4 8], ...
                'display', dispRate);
            n1 = connectAdd(s, n1, 1);

            if oto
                n2 = connectAdd(n1, n2, [], 'weight', w, 'type', connectionTypes{ct}, ...
                    'learn', learningParams.lambda(lr), learningParams.mu1(lr),...
                    learningParams.mu2(lr), learningParams.ceps(lr),...
                    learningParams.kappa(lr), ...
                    'display', dispRate, 'phasedisp', 'save', 1000);
            else
                n2 = connectAdd(n1, n2, [], 'weight', w, 'type', connectionTypes{ct}, 'no11', ...
                    'learn', learningParams.lambda(lr), learningParams.mu1(lr),...
                    learningParams.mu2(lr), learningParams.ceps(lr),...
                    learningParams.kappa(lr), ...
                    'display', dispRate, 'phasedisp', 'save', 1000);
            end
            
            M = modelMake(@zdot, @cdot, s, n1, n2);
            
            %% Run the network
            tic
            M = odeRK4fs(M);
            runTimes(ct,lr) = toc;
            Z = M.n{1}.Z;
            
            if oto
                nansPresent1to1(ct,lr) = any(isnan(Z(:)));
                if (any(isnan(Z(:))))
                    disp('Warning, NaNs present');
                else
                    disp('OK, no NaNs');
                end
            else
                nansPresentno11(ct,lr) = any(isnan(Z(:)));
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

if sum(nansPresentno11(:))>0 || sum(nansPresent1to1(:))>0
    disp(' ')
    disp('NaNs present - check above output results')
else
    disp(' ')
    disp('No NaNs present in any of the above simulations')
end