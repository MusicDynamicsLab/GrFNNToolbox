%% specAndAutocorr
%  specAndAutocorr(y,NFFT,fs,varargin)
%
%  Plotting function. Produces a figure with a spectrogram and autocorrelogram
%  of signal vector y. Input arguments for fft length NFFT and sampling 
%  frequency fs are required. Optional inputs must come in order afterwards, 
%  including portion, N, and windowstep. portion is a two element vector of 
%  percentages indicating the range of frequencies, from 0 Hz to the Nyquist, 
%  that should be displayed. N is the length of each window in samples. 
%  windowStep is the length by which to step to the next window, in samples. 

function specAndAutocorr(y,NFFT,fs,varargin)

if nargin<3,error('specAndAutocorr needs at least 3 inputs: signal vector, fft size, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
end
y=y(:);
if isempty(varargin)
    portion=[0 100];
    N=ceil(length(y)^(1/1.4));
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==1
    portion=varargin{1};
    N=ceil(length(y)^(1/1.4));
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==2
    portion=varargin{1};
    N=varargin{2};
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==3
    portion=varargin{1};
    N=varargin{2};
    windowStep=varargin{3};
else
    error('specAndAutocorr only takes 6 inputs')
end
if ~isreal(y)
    warning('Input signal is complex; only real portion taken')
    y=real(y);
end

t=(0:length(y))/fs;

f = fs/2*linspace(0,1,NFFT/2+1);

f1 = floor((portion(1)/100)*length(f))+1;
f2 = floor((portion(2)/100)*length(f));
fRange = [f(f1) f(f2)];



screen = get(0,'screensize');
screen(1:2)=screen(3:4)/4;
screen(3:4)=screen(3:4)/2;
figure;set(gcf,'position',screen);
setappdata(gcf, 'SubplotDefaultAxesLocation', [.08 .085 .91 .87]);

subplot(1,2,1)

struct = mdlSpec(y,NFFT,fs,portion,N,windowStep);
Sfft = struct.spec;

i1 = floor((portion(1)/200)*size(Sfft,1))+1;
i2 = floor((portion(2)/200)*size(Sfft,1));
imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))+1));

title('Spectrogram');
xlabel('Time (sec)');ylabel('Frequency (Hz)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');

subplot(1,2,2)

mdlAutocorr(y,fs,N,windowStep);
title('Autocorrelogram');
xlabel('Time (sec)');ylabel('Lags (ms)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Correlation (-1 to 1)');
drawnow;