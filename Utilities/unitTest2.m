%% unitTest2.m
%
% A simple afferent chain network with fixed connections
%

connectionTypes = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};

aLin  = -1; aCrit  = 0; aCritDetune =  0; aLC = 1;   aDLC = -1;
b1Lin =  0; b1Crit =-1; b1CritDetune= -1; b1LC=-1;   b1DLC = 3;
b2Lin =  0; b2Crit =-1; b2CritDetune= -1; b2LC=-1e3; b2DLC =-1;
d1Lin =  0; d1Crit = 0; d1CritDetune=  1; d1LC= 0;   d1DLC = 0;
d2 = 0;
eLin  =  1; eCrit  = 1; eCritDetune =  1; eLC = 1;   eDLC  = 1;

params = struct('alpha', [aLin, aCrit, aCritDetune, aLC, aDLC],...
    'beta1', [b1Lin, b1Crit, b1CritDetune, b1LC, b1DLC],...
    'beta2', [b2Lin, b2Crit, b2CritDetune, b2LC, b2DLC],...
    'delta1', [d1Lin, d1Crit, d1CritDetune, d1LC, d1DLC],...
    'delta2', [d2, d2, d2, d2, d2],...
    'eps', [eLin, eCrit, eCritDetune, eLC, eDLC]);

Fs = 160;
w = 1;
dispRate = 100;

%% Make the model
for ct = 1:length(connectionTypes)
    disp(' ')
    disp(['---Connection type: ' connectionTypes{ct} '---']);
    for oto = 1:-1:0 % run with and without one to one connections
        if oto==1
            disp(' ')
            disp('including 1:1 connections')
        elseif ~oto && (ct==1 || ct==3)
            disp(' ')
            disp('testing no 1:1 not necessary for this connection type')
            continue;
        else
            disp(' ')
            disp('no 1:1 connections')
        end
        for pr = 1:5
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
            
            Nosc = 201;
            if ct==4 % if 3freqall
                Nosc = (Nosc-1)/4 + 1;
            end
            
            s = stimulusMake(1, 'fcn', [0 10], Fs, {'exp'}, [1 4], .25, 0, 'ramp', 0.01, 1, ...
                'display', dispRate);
            n1 = networkMake(1, 'hopf', params.alpha(2), params.beta1(2), params.beta2(2),...
                params.delta1(2), params.delta2(2), params.eps(2),...
                'log', .5, 8, Nosc, 'save', 1, 'Tick', [.5 1 2 4 8],...
                'display', dispRate);
            
            n2 = networkMake(2, 'hopf', params.alpha(pr), params.beta1(pr), params.beta2(pr),...
                params.delta1(pr), params.delta2(pr), params.eps(pr),...
                'log', .5, 8, (Nosc-1)/2+1, 'save', 1, 'Tick', [.5 1 2 4 8],...
                'display', dispRate);
            
            if ct==3 || ct==4 % if 3freq or 3freqall
                C = .005;
            else
                C = connectMake(n1, n2, 'full', 0.005);
            end
            
            n1 = connectAdd(s, n1, 1); % '1freq' connection type by default
            
            if oto
                n2 = connectAdd(n1, n2,  C, 'weight', w, 'type', connectionTypes{ct});
            else
                n2 = connectAdd(n1, n2,  C, 'weight', w, 'type', connectionTypes{ct}, 'no11');
            end
            
            M = modelMake(@zdot, @cdot, s, n1, n2);
            
            %% Run the network
            tic
            M = odeRK4fs(M);
            runTimes(ct,pr) = toc;
            
            Z = [M.n{1}.Z; M.n{2}.Z];
            
            if oto
                nansPresent1to1(ct,pr) = any(isnan(Z(:)));
                if (any(isnan(Z(:))))
                    disp('Warning, NaNs present');
                else
                    disp('OK, no NaNs');
                end
            else
                nansPresentno11(ct,pr) = any(isnan(Z(:)));
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