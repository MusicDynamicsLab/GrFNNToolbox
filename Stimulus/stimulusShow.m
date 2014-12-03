%% stimulusShow
%  1 - plot
%  2 - spectrogram (requires f, tck, and tckl variables)
%
%  stimulus channels among rows, samples along columns     
%  ie, for a two-channel input with 40,000 samples, 
%  size(s.x) should be [2 40000]

%%
function stimulusShow(s, windows, f, tck, tckl)

if nargin < 2
    windows = [1];
end

if ismember(1,windows) & isfield(s,'x')

    figure(01)

    for c = 1:size(s.x,1)
        subplot(size(s.x,1),1,c)
        set(gca, 'XLim', s.ts([1 end]));
        plot(s.t, real(s.x(c,:)), 'k')
        %title(['Channnel ', num2str(c)])
        set(gca, 'XLim', [s.t(1) s.t(end)])
        set(gca, 'YLim', [-1 1] * 1.2*(max(abs(s.x(:)))+eps))
        
        xlabel('Time')
        ylabel('Amplitude')
        title(['Stimulus ', 'Channel ', num2str(c)]);
        grid on
        zoom xon
        
    end
end

if ismember(2,windows) & nargin > 2

    figure(02)

    Z = spectrogram(real(s.x),8*s.fs,round(0.95*8*s.fs), f, s.fs);
    imagesc(s.t, 1:length(f), abs(Z)); axis xy
    set(gca, 'XLim', s.ts([1 end]));
    set(gca, 'YTick', tck, 'YTickLabel', tckl)
    colormap(flipud(hot))
    colorbar
    title('Spectrogram of Stimulus');
    grid on
    zoom xon

end
