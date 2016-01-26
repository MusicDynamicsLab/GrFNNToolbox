%% odeRK4fs
%   M = odeRK4fs(M)
%
%  Fixed-step 4th-order Runge-Kutta ODE numerical integration.
%  Steps by *direct indexing* of stimulus vector.
%
%  Params
%   Model
%
%  Output
%   M - Model

%%
function M = odeRK4fs(M, varargin)

load('MyColormaps', 'IF_colormap');
circular = IF_colormap;

zfun = M.dotfunc;
cfun = M.cfun;
ispan = [1 length(M.t)];

%% Error checking
if( ~isa(zfun,'function_handle') )
    error('odeRK4fs: odefun param must be a function handle');
end

step = single(1);
h = single(M.dt);                   % step size
stimList = M.stimList;
netList = M.netList;

%% Display stimulus and initial conditions if dStep > 0
for sx = stimList
    if M.s{sx}.dStep
        M.s{sx}.bH = stimulusDisplay(M.s{sx}, 0, M.t(1));
    end
end

for nx = netList
    if M.n{nx}.dStep
        [M.n{nx}.nH, M.n{nx}.tH] = networkDisplay(M.n{nx}, 0, M.t(1));
    end
    for cx = M.n{nx}.conLearn
        if M.n{nx}.con{cx}.dStep
            [M.n{nx}.con{cx}.aH, M.n{nx}.con{cx}.atH, M.n{nx}.con{cx}.pH, M.n{nx}.con{cx}.ptH] ...
                = connectionDisplay(M.n{nx}.con{cx}, 0, M.t(1), circular);
        end
    end
end

%% Integration loop
for ix = ispan(1) : step : ispan(2)-step
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
            M.n{nx}.k{kx} = h*zfun(M, nx);

            %% ... and for each learned connection to the network
            for cx = M.n{nx}.conLearn
                M.n{nx}.con{cx}.k{kx} = h*cfun(M, nx, cx);
            end
        end
        
        %% Update z, C and ind for the next k-step
        switch kx
            case 1
                for nx = netList
                    M.n{nx}.zPrev = M.n{nx}.z;
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{1}/2;
                    for cx = M.n{nx}.conLearn
                        M.n{nx}.con{cx}.CPrev = M.n{nx}.con{cx}.C;
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{1}/2;
                    end
                end
                ind = ix + step/2; % time step for k2 and k3
            case 2
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{2}/2;
                    for cx = M.n{nx}.conLearn
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{2}/2;
                    end
                end
            case 3
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{3};
                    for cx = M.n{nx}.conLearn
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{3};
                    end
                end
                ind = ix + step; % time step for k4
            case 4
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + ...
                        (M.n{nx}.k{1} + 2*M.n{nx}.k{2} + 2*M.n{nx}.k{3} + M.n{nx}.k{4})/6;
                    if M.n{nx}.sStep && ~mod(ix, M.n{nx}.sStep)
                        M.n{nx}.Z(:,ix/M.n{nx}.sStep+1) = M.n{nx}.z;
                    end
                    if M.n{nx}.dStep && ~mod(ix, M.n{nx}.dStep)
                        networkDisplay(M.n{nx}, ix, M.t(ix));
                    end
                    for cx = M.n{nx}.conLearn
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + ...
                            (M.n{nx}.con{cx}.k{1} + 2*M.n{nx}.con{cx}.k{2} + 2*M.n{nx}.con{cx}.k{3} + M.n{nx}.con{cx}.k{4})/6;
                        if M.n{nx}.con{cx}.sStep && ~mod(ix, M.n{nx}.con{cx}.sStep)
                            M.n{nx}.con{cx}.C3(:,:,ix/M.n{nx}.con{cx}.sStep+1) = M.n{nx}.con{cx}.C;
                        end
                        if M.n{nx}.con{cx}.dStep && ~mod(ix, M.n{nx}.con{cx}.dStep)
                            connectionDisplay(M.n{nx}.con{cx}, ix, M.t(ix), circular);
                        end
                    end
                end
                
                for sx = stimList
                    if M.s{sx}.dStep && ~mod(ix, M.s{sx}.dStep)
                        stimulusDisplay(M.s{sx}, ix, M.t(ix));
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


