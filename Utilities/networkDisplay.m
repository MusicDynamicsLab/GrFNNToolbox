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