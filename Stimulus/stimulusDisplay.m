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