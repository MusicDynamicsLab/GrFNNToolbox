%% function: Displays instantaneous network state
function networkLiveDisplay(M, nx, ix, t)

persistent networkDispMap;
net = M.n{nx};

if ix == 0    
    if isfield(net,'nAx') && ishghandle(net.nAx)
        networkData.nAx = net.nAx;
        axes(net.nAx)
    else
        figure(10000+1000*net.id)
	end
    
	networkData.nH = plot(net.f, abs(net.z), '.-');  % nH: lineseries object handle
    switch net.nFspac
        case 1 % linear spacing
            networkData.nH = plot(net.f, abs(net.z), '.-');  % nH: lineseries object handle
        case 2 % log spacing
            networkData.nH = semilogx(net.f, abs(net.z), '.-');  % nH: lineseries object handle
    end
    
    if isfield(net,'title') && ischar(net.title)
        title(net.title)
        networkData.tH = [];
    else
        networkData.tH = title(sprintf('Amplitudes of oscillators in network %d (t = %.1fs)', net.id, t));
    end
    
    xlabel('Oscillator natural frequency (Hz)')
    ylabel('Amplitude')
    set(gca, 'XLim',[min(net.f) max(net.f)])
    set(gca, 'YLim', [0 1/sqrt(net.e)])
    if isfield(net,'tick') && ~isempty(net.tick)
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