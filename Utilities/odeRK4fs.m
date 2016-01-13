%% odeRK4fs
%   M = odeRK4fs(M,s)
%
%  Fixed-step 4th-order Runge-Kutta ODE numerical integration.
%  Steps by *direct indexing* of stimulus vector.
%
%  Params
%   Model, stimulus
%
%  Output
%   M - Model

%%
function M = odeRK4fs(M, varargin)

% clear global M
% global M circular

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
h = single(M.dt);                   % For variable step size, else h = dt;
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

        %% ... for each network
        for nx = netList
            M.n{nx}.k{kx} = h*zfun(M, ind, nx);

            %% ... and for each learned connection to the network
            for cx = M.n{nx}.conLearn
                M.n{nx}.con{cx}.k{kx} = h*cfun(M, ind, nx, cx);
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

%% function: Displays stimulus and progress bar
function bH = stimulusDisplay(stim, ix, t)
% global M 
if ix == 0
    sx = stim.id;
    if isfield(stim,'sAx') && ishghandle(stim.sAx)
        axes(stim.sAx)
    else
        figure(10000+sx)
    end
    
    plot(stim.t, real(stim.x(stim.dispChan,:)), 'k')
    hold on
    ylimit = get(gca,'YLim');
    bH = plot([1 1]*t, ylimit, 'r'); % handle for progress bar
    set(gca,'YLim',ylimit) % need to do this for lastest matlab
    set(gca,'XLim',[min(stim.t) max(stim.t)])
    hold off
    if isfield(stim,'title') && ischar(stim.title)
        title(stim.title)
    else
        title(['Stimulus ' num2str(sx) ', Channel ', num2str(stim.dispChan)])
    end
    xlabel('Time')
    ylabel('Real part')
    grid
    
else
    set(stim.bH, 'XData', [1 1]*t)
end

drawnow

end

%% function: Displays instantaneous network state
function [nH, tH] = networkDisplay(net, ix, t)
% global M
nx = net.id;
if ix == 0
    if isfield(net,'nAx') && ishghandle(net.nAx)
        axes(net.nAx)
    else
        figure(10000+1000*nx)
    end
    
    switch net.fspac
        case 'log'
            nH = semilogx(net.f, abs(net.z), '.-');  % nH: lineseries object handle
        case 'lin'
            nH = plot(net.f, abs(net.z), '.-');  % nH: lineseries object handle
    end
    if isfield(net,'title') && ischar(net.title)
        title(net.title)
        tH = [];
    else
        tH = title(sprintf('Amplitudes of oscillators in network %d (t = %.1fs)', nx, t));
    end
    xlabel('Oscillator natural frequency (Hz)')
    ylabel('Amplitude')
    set(gca, 'XLim',[min(net.f) max(net.f)])
    set(gca, 'YLim', [0 1/sqrt(net.e)])
    if ~isempty(net.tick)
        set(gca, 'XTick', net.tick)
    end
    
    grid
else
    set(net.nH, 'YData', abs(net.z))
    if ~isempty(net.tH)
        set(net.tH, 'String', sprintf('Amplitudes of oscillators in network %d (t = %.1fs)', nx, t))
    end
end

drawnow

end

%% function: Displays instantaneous connection state
function [aH, atH, pH, ptH] = connectionDisplay(con, ix, t, cmap)
% global M circular
nx = con.target;
cx = con.id;

if ix == 0
    if isfield(con,'aAx') && ishghandle(con.aAx)
        axes(con.aAx)
    elseif ~ishghandle(10000+1000*nx+100*cx)
        figure(10000+1000*nx+100*cx)
        set(gcf, 'Position', [2 550 500 400])
    else
        figure(10000+1000*nx+100*cx)
    end
    is3freq = con.nType == 3 || con.nType == 4;   % if 3freq or 3freqall
    
    if size(con.C, 2) == 1 % if connection is a column vector
        aH = plot(con.targetAxis, abs(con.C), '.-');
        xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.target))
        ylabel('Connection amplitude')
        set(gca, 'XLim',[min(con.targetAxis) max(con.targetAxis)])
        if abs(con.e)
            set(gca, 'YLim', [0 1/sqrt(con.e)])
        end
        if strcmp(con.targetAxisScale, 'log')
            set(gca, 'XScale', 'log')
        end
        if ~isempty(con.targetAxisTick)
            set(gca, 'XTick', con.targetAxisTick);
        end
        
    else    % if connection is a matrix
        if is3freq || con.nSourceClass == 1 % 3freq or stimulus source
            f1 = [1 size(con.C, 2)];
        else
            f1 = getLim(con.sourceAxis, con.sourceAxisScale);
        end
        f2 = getLim(con.targetAxis, con.targetAxisScale);
        aH = imagesc(f1, f2, abs(con.C));
        colormap(flipud(hot)); colorbar;
        ylabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.target))
        if abs(con.e)
            set(gca, 'CLim', [.001 1/sqrt(con.e)])
        end
        if strcmp(con.sourceAxisScale, 'log') && ~is3freq
            set(gca, 'XScale', 'log')
        end
        if strcmp(con.targetAxisScale, 'log')
            set(gca, 'YScale', 'log')
        end
        if ~isempty(con.sourceAxisTick) && ~is3freq
            set(gca, 'XTick', con.sourceAxisTick);
        end
        if ~isempty(con.targetAxisTick)
            set(gca, 'YTick', con.targetAxisTick);
        end
        if is3freq
            xlabel(sprintf('Oscillator pair: Network %d', con.source))
        elseif con.nSourceClass == 1 && size(con.C, 2) > 1  % channelized stimulus
            xlabel(sprintf('Channel: Stimulus %d', con.source))
        else
            xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.source))
        end
    end
    
    if isfield(con,'titleA')
        title(con.titleA)
        atH = [];
    else
        atH = title(sprintf('Amplitudes of connection matrix %d to network %d (t = %.1fs)', cx, nx, t));
    end
    grid on
    
    if con.phaseDisp
        if isfield(con,'pAx') && ishghandle(con.pAx)
            axes(con.pAx)
        elseif ~ishghandle(10000+1000*nx+100*cx+1)
            figure(10000+1000*nx+100*cx+1);
            set(gcf, 'Position', [500 550 500 400])
        else
            figure(10000+1000*nx+100*cx+1)
        end
        
        if size(con.C, 2) == 1
            pH = plot(con.targetAxis, angle(con.C), '.');
            xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.target))
            ylabel('Connection phase')
            set(gca, 'XLim',[min(con.targetAxis) max(con.targetAxis)])
            set(gca, 'YLim', [-pi pi])
            set(gca, 'YTick', [-pi, -pi/2, 0, pi/2, pi])
            set(gca, 'YTickLabel', {'-pi  ', '-pi/2', ' 0  ', ' pi/2', ' pi  '})
            if strcmp(con.targetAxisScale, 'log')
                set(gca, 'XScale', 'log')
            end
            if ~isempty(con.targetAxisTick)
                set(gca, 'XTick', con.targetAxisTick);
            end

        else
            pH = imagesc(f1, f2, angle(con.C));
            colormap(gca, cmap);
            cb = colorbar;
            ylabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.target))
            set(cb, 'YTick',      [-pi, -pi/2, 0, pi/2, pi])
            set(cb, 'YTickLabel', {'-pi  ', '-pi/2', ' 0  ', ' pi/2', ' pi  '})
            set(gca, 'CLim', [-pi pi])
            if strcmp(con.sourceAxisScale, 'log') && ~is3freq
                set(gca, 'XScale', 'log')
            end
            if strcmp(con.targetAxisScale, 'log')
                set(gca, 'YScale', 'log')
            end
            if ~isempty(con.sourceAxisTick) && ~is3freq
                set(gca, 'XTick', con.sourceAxisTick);
            end
            if ~isempty(con.targetAxisTick)
                set(gca, 'YTick', con.targetAxisTick);
            end
            if is3freq
                xlabel(sprintf('Oscillator pair: Network %d', con.source))
            elseif con.nSourceClass == 1 && size(con.C, 2) > 1  % channelized stimulus
                xlabel(sprintf('Channel: Stimulus %d', con.source))
            else
                xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.source))
            end
        end
        
        if isfield(con,'titleP') && ischar(con.titleP)
            title(con.titleP)
            ptH = [];
        else
            ptH = title(sprintf('Phases of connection matrix %d to network %d (t = %.1fs)', cx, nx, t));
        end
        grid on
        
    else    % if no phase display
        pH = [];
        ptH = [];
    end

else    % nonzero ix
    if size(con.C, 2) == 1
        set(con.aH, 'YData', (abs(con.C)))
    else
        set(con.aH, 'CData', (abs(con.C)))
    end
    if ~isempty(con.atH)
        set(con.atH, 'String', sprintf('Amplitudes of connection matrix %d to network %d (t = %.1fs)', cx, nx, t))
    end
    if con.phaseDisp
        if size(con.C, 2) == 1
            set(con.pH, 'YData', angle(con.C))
        else
            set(con.pH, 'CData', angle(con.C))
        end
        if ~isempty(con.ptH)
            set(con.ptH, 'String', sprintf('Phases of connection matrix %d to network %d (t = %.1fs)', cx, nx, t))
        end
    end
end

drawnow

end
