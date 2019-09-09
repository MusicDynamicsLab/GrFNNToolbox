%% odeRK4fsNorm
%   M = odeRK4fsNorm(M)
%
%  Integrates a model, M, using fixed-step 4th-order Runge-Kutta method.
%  Steps by direct indexing of stimulus vector.
%  Works with zfun that passes a netwokrk object (e.g. zdotNorm.m)
%

%%
function M = odeRK4fsNorm(M)

zfun = M.zfun;
cfun = M.cfun;
iSpan = M.iSpan;
iStep = M.iStep;
h = M.dt;                   % step size
stimList = M.stimList;
netList = M.netList;

%% Display stimulus and initial conditions if dStep > 0
for sx = stimList
    if M.s{sx}.dStep
        stimulusLiveDisplay(M, sx, 0, M.t(1));
    end
end

for nx = netList
    if M.n{nx}.dStep
        networkLiveDisplay(M, nx, 0, M.t(1));
    end
    for cx = M.n{nx}.learnList
        if M.n{nx}.con{cx}.dStep
            connectionLiveDisplay(M, nx, cx, 0, M.t(1));
        end
    end
end

%% Integration loop
for ix = iSpan(1) : iStep : iSpan(2)-iStep
    ind = ix; % time step for which to calculate k1
    
    %% Get Runge-Kutta k-values
    for kx = 1:4
        
        %% First update stimulus values
        if kx == 2 || kx == 4
            for sx = stimList
                M.s{sx}.z = stimulusRun(M.s{sx}, ind);
            end
        end

        %% Get k-values for each network
        for nx = netList
%             M.n{nx}.k{kx} = h*zfun(M, nx);
            M.n{nx} = zfun(M, nx, kx); % zfun now passes a network

            %% ... and for each learned connection to the network
            for cx = M.n{nx}.learnList
                M.n{nx}.con{cx}.k{kx} = h*cfun(M, nx, cx);
            end
        end
        
        %% Update z, C and ind for the next k-step
        switch kx
            case 1
                for nx = netList
                    M.n{nx}.zPrev = M.n{nx}.z;
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{1}/2;
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.CPrev = M.n{nx}.con{cx}.C;
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{1}/2;
                    end
                end
                ind = ix + iStep/2; % time step for k2 and k3
            case 2
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{2}/2;
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{2}/2;
                    end
                end
            case 3
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{3};
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{3};
                    end
                end
                ind = ix + iStep; % time step for k4
            case 4
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + ...
                        (M.n{nx}.k{1} + 2*M.n{nx}.k{2} + 2*M.n{nx}.k{3} + M.n{nx}.k{4})/6;
                    if M.n{nx}.sStep && ~mod(ix, M.n{nx}.sStep)
                        M.n{nx}.Z(:,ix/M.n{nx}.sStep+1) = M.n{nx}.z;
                    end
                    if M.n{nx}.dStep && ~mod(ix, M.n{nx}.dStep)
                        networkLiveDisplay(M, nx, ix, M.t(ix));
                    end
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + ...
                            (M.n{nx}.con{cx}.k{1} + 2*M.n{nx}.con{cx}.k{2} + 2*M.n{nx}.con{cx}.k{3} + M.n{nx}.con{cx}.k{4})/6;
                        if M.n{nx}.con{cx}.sStep && ~mod(ix, M.n{nx}.con{cx}.sStep)
                            M.n{nx}.con{cx}.C3(:,:,ix/M.n{nx}.con{cx}.sStep+1) = M.n{nx}.con{cx}.C;
                        end
                        if M.n{nx}.con{cx}.dStep && ~mod(ix, M.n{nx}.con{cx}.dStep)
                            connectionLiveDisplay(M, nx, cx, ix, M.t(ix));
                        end
                    end
                end
                
                for sx = stimList
                    if M.s{sx}.dStep && ~mod(ix, M.s{sx}.dStep)
                        stimulusLiveDisplay(M, sx, ix, M.t(ix));
                    end
                end
        end
    end
end

%% function: computes one 4th-order Runge-Kutta step
%   See http://www.physics.utah.edu/~detar/phys6720/handouts/ode/ode/node6.html
%
%  k1 = h*f(ti, yi)
%  k2 = h*f(ti+h/2, yi+k1/2)
%  k3 = h*f(ti+h/2, yi+k2/2)
%  k4 = h*f(t(i+1), yi+k3)
%
%  y(i+1) = yi + 1/6*(k1 + 2*k2 + 2*k3 + k4)


