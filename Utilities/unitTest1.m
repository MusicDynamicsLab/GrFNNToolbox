%% unitTest1.m
%
% A one layer network driven with a sinusoidal input. Several parameter
% sets are provided for experimentation with different types of intrinsic
% oscillator dynamics.

connectionTypes = {'1freq', 'all2freq', 'allfreq', 'active'};

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

fs = 160;
dispRate = 10;

numConnectionTypes = length(connectionTypes);
count = 0;

for ct = 1:numConnectionTypes      % loop over connection types
    disp(' ')
    disp(['---Connection type: ' connectionTypes{ct} '---']);
    for oto = 1:-1:0 % run with and without one to one connections
        if oto==1
            disp(' ')
            disp('including 1to1 connections')
        elseif ~oto && ct==1
            disp(' ')
            disp('testing no11 not necessary for this connection type')
            continue;
        else
            disp(' ')
            disp('no 1to1 connections')
        end
        for pr = 1:5 % parameter regime loop
            count = count+1;
            switch pr
                case 1
                    fprintf('- Linear regime running...')
                case 2
                    fprintf('- Critical regime running...')
                case 3
                    fprintf('- Critical regime with detuning running...')
                case 4
                    fprintf('- Limit cycle regime running...')
                case 5
                    fprintf('- Double limit cycle regime running...')
            end
            %% Make the model
            s = stimulusMake(1, 'fcn', [0 50], fs, {'exp'}, [2], .25, 0, ...
                'ramp', 0.02, 1, 'display', dispRate, 'InputType', connectionTypes{ct});
                        
            n = networkMake(1, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 8, 201, 'save', 1, ...
                'display', dispRate, 'Tick', [.5 1 2 4 8]);
            
            C = ones(n.N, s.N);
            
            if oto
                n = connectAdd(s, n, C, 'type', connectionTypes{ct});
            elseif ~oto && ct~=1
                n = connectAdd(s, n, C, 'type', connectionTypes{ct}, 'no11');
            end
            
            %
            % Run the integrator
            %
            
            M = modelMake(@zdot, @cdot, s, n);
            
            tic
            M = odeRK4fs(M);
            runTimes(count) = toc;
            
            Z = M.n{1}.Z;
            
            
            if oto
                nansPresent1to1(pr,ct) = any(isnan(Z(:)));
                if (any(isnan(Z(:))))
                    disp('Warning, NaNs present');
                else
                    disp('OK, no NaNs');
                end
            else
                nansPresentno11(pr,ct) = any(isnan(Z(:)));
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
