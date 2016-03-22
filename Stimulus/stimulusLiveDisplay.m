%% function: Displays stimulus and progress bar
function stimulusLiveDisplay(M, sx, ix, t)

persistent stimulusDispMap;
stim = M.s{sx};

if ix == 0
    
    sx = stim.id;
    if isfield(stim,'sAx') && ishghandle(stim.sAx)
        stimulusData.sAx = stim.sAx;
        axes(stim.sAx)
    else
        figure(10000+sx)
    end
    
    plot(stim.t, real(stim.x(stim.dispChan,:)), 'k')
    hold on
    ylimit = get(gca,'YLim');
    stimulusData.bH = plot([1 1]*t, ylimit, 'r'); % handle for progress bar
    set(gca,'YLim',ylimit) % need to do this for latest matlab
    set(gca,'XLim',[min(stim.t) max(stim.t)])
    hold off
    if isfield(stim,'title') && ischar(stim.title)
        stimulusData.title = title(stim.title);
    else
        stimulusData.title = ...
            title(['Stimulus ' num2str(sx) ', Channel ', num2str(stim.dispChan)]);
    end
    xlabel('Time')
    ylabel('Real part')
    grid
    stimulusDispMap{stim.id} = stimulusData;
else
    stimulusData = stimulusDispMap{stim.id};
    set(stimulusData.bH, 'XData', [1 1]*t)
end

drawnow