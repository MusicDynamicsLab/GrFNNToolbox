%% function: Displays instantaneous network state
function [nH, tH] = networkDisplay(net, ix, t)

persistent networkDispMap;

if ix == 0
    if (size(networkDispMap,2) < net.id)
        networkData.id = net.id;
        networkDispMap{net.id} = networkData;
    end
    
    if isfield(net,'nAx') && ishghandle(net.nAx)
        networkData.nAx = axes(net.nAx);
    else
        figure(10000+1000*net.id)
    end
    
    switch net.fspac
        case 1 % linear spacing
            networkData.nH = semilogx(net.f, abs(net.z), '.-');  % nH: lineseries object handle
            nH = networkData.nH;
        case 2 % log spacing
            networkData.nH = plot(net.f, abs(net.z), '.-');  % nH: lineseries object handle
            nH = networkData.nH;
    end
    if isfield(net,'title') && ischar(net.title)
        title(net.title)
        networkData.tH = [];
        tH = networkData.tH;
    else
        networkData.tH = title(sprintf('Amplitudes of oscillators in network %d (t = %.1fs)', net.id, t));
        tH = networkData.tH;
    end
    xlabel('Oscillator natural frequency (Hz)')
    ylabel('Amplitude')
    set(gca, 'XLim',[min(net.f) max(net.f)])
    set(gca, 'YLim', [0 1/sqrt(net.e)])
    if ~isempty(net.tick)
        set(gca, 'XTick', net.tick)
    end
    
    grid
    networkDispMap{net.id} = networkData;
else
    networkData = networkDispMap{net.id};
    
    set(networkData.nH, 'YData', abs(net.z))
    if ~isempty(networkData.tH)
        set(networkData.tH, 'String',...
            sprintf('Amplitudes of oscillators in network %d (t = %.1fs)', net.id, t))
    end
end

drawnow