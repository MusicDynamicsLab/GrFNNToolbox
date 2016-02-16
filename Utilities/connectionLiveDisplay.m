%% function: Displays instantaneous connection state
function connectionLiveDisplay(con, ix, t, cmap)

nx = con.target;
cx = con.id;
persistent connectionDispMap;

if ix == 0
    if isfield(con,'aAx') && ishghandle(con.aAx)
        connectionData.aAx = axes(con.aAx);
    elseif ~ishghandle(10000+1000*nx+100*cx)
        figure(10000+1000*nx+100*cx)
        set(gcf, 'Position', [2 550 500 400])
    else
        figure(10000+1000*nx+100*cx)
    end
    is3freq = con.nType == 3 || con.nType == 4;   % if 3freq or 3freqall
    
    if con.sourceN == 1 % if connection is a column vector
        connectionData.aH = plot(con.targetAxis, abs(con.C), '.-');
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
            f1 = [1 con.sourceN];
        else
            f1 = getLim(con.sourceAxis, con.nSourceAxisScale);
        end
        f2 = getLim(con.targetAxis, con.nTargetAxisScale);
        connectionData.aH = imagesc(f1, f2, abs(con.C));
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
            if con.nSourceClass == 1 % stimulus source
                xlabel(sprintf('Stimulus pair: Stimulus %d', con.source))
            else
                xlabel(sprintf('Oscillator pair: Network %d', con.source))
            end
        elseif con.nSourceClass == 1 && con.sourceN > 1  % channelized stimulus
            xlabel(sprintf('Channel: Stimulus %d', con.source))
        else
            xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.source))
        end
    end
    
    if isfield(con,'titleA')
        title(con.titleA)
        connectionData.atH = [];
    else
        connectionData.atH = title(sprintf('Amplitudes of connection matrix %d to network %d (t = %.1fs)', cx, nx, t));
    end
    grid on
    
    if con.phaseDisp
        if isfield(con,'pAx') && ishghandle(con.pAx)
            connectionData.pAx = axes(con.pAx);
        elseif ~ishghandle(10000+1000*nx+100*cx+1)
            figure(10000+1000*nx+100*cx+1);
            set(gcf, 'Position', [500 550 500 400])
        else
            figure(10000+1000*nx+100*cx+1)
        end
        
        if con.sourceN == 1
            connectionData.pH = plot(con.targetAxis, angle(con.C), '.');
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
            connectionData.pH = imagesc(f1, f2, angle(con.C));
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
                if con.nSourceClass == 1 % stimulus source
                    xlabel(sprintf('Stimulus pair: Stimulus %d', con.source))
                else
                    xlabel(sprintf('Oscillator pair: Network %d', con.source))
                end
            elseif con.nSourceClass == 1 && con.sourceN > 1  % channelized stimulus
                xlabel(sprintf('Channel: Stimulus %d', con.source))
            else
                xlabel(sprintf('Oscillator natural frequency (Hz): Network %d', con.source))
            end
        end
        
        if isfield(con,'titleP') && ischar(con.titleP)
            title(con.titleP)
            connectionData.ptH = [];
        else
            connectionData.ptH = title(sprintf('Phases of connection matrix %d to network %d (t = %.1fs)', cx, nx, t));
        end
        grid on
        
    else    % if no phase display
        connectionData.pH = [];
        connectionData.ptH = [];
    end
    connectionDispMap{con.id} = connectionData;

else    % nonzero ix
    connectionData = connectionDispMap{con.id};
    if con.sourceN == 1
        set(connectionData.aH, 'YData', (abs(con.C)))
    else
        set(connectionData.aH, 'CData', (abs(con.C)))
    end
    if ~isempty(connectionData.atH)
        set(connectionData.atH, 'String', sprintf('Amplitudes of connection matrix %d to network %d (t = %.1fs)', cx, nx, t))
    end
    if con.phaseDisp
        if con.sourceN == 1
            set(connectionData.pH, 'YData', angle(con.C))
        else
            set(connectionData.pH, 'CData', angle(con.C))
        end
        if ~isempty(connectionData.ptH)
            set(connectionData.ptH, 'String', sprintf('Phases of connection matrix %d to network %d (t = %.1fs)', cx, nx, t))
        end
    end
end

drawnow