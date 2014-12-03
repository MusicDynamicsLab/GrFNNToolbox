function allMyFreqs(y,NFFT,fs,varargin)

if nargin<3,error('myFreqs needs at least 3 inputs: signal vector, fft size, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
end
y=y(:);
loop=0;
dbg=0;
if isempty(varargin)
    portion=[0 100];
    N=floor(fs/20);
    windowStep=50;
elseif length(varargin)==1
    portion=varargin{1};
    N=floor(fs/20);
    windowStep=50;
elseif length(varargin)==2
    portion=varargin{1};
    N=varargin{2};
    windowStep=50;
elseif length(varargin)==3
    portion=varargin{1};
    N=varargin{2};
    windowStep=varargin{3};
elseif length(varargin)==4 && strcmpi(varargin{4},'loop')
    portion=varargin{1};
    N=varargin{2};
    windowStep=varargin{3};
    loop=1;
elseif length(varargin)==5 && strcmpi(varargin{5},'dbg')
    portion=varargin{1};
    N=varargin{2};
    windowStep=varargin{3};
    if strcmpi(varargin{4},'loop'), loop=1;end
    dbg=1;
else
    error('myFreqs only takes 6 inputs')
end
if ~isreal(y)
    warning('Input signal is complex; only real portion taken')
    y=real(y);
end

if dbg, tic; end;

t=(0:length(y))/fs;

f = fs/2*linspace(0,1,NFFT/2+1);
% if mod(N,2),N=N+1;end
n = 0:N-1;
window = exp(-.5*((n-(N-1)/2)/(.4*(N-1)/2)).^2)';
% window = gausswin(N)';
% window = ones(N,1);
% window = hamming(N)';
% figure;plot(n,window);drawnow;

% x = zeros(length(y)+2*N,1);
% x(N:length(y)+N-1) = y;

if loop
    
    %     Tstep = floor((length(y)+N)/windowStep);
    Tstep = floor((length(y)-N)/windowStep);
    Sfft = zeros(length(f),Tstep);
    
    count = 1;
    
    %     for i = 1:windowStep:length(y)+N
    for i = 1:windowStep:length(y)-N
        
        %         Y = fft(window.*x(i:i+N-1),NFFT);
        Y = fft(window.*y(i:i+N-1),NFFT);
        Sfft(:,count) = Y(1:NFFT/2+1).';
        
        count = count + 1;
    end;if dbg, toc; end;
    
else
    
    %     Sfft = fft(repmat(window,1,ceil((length(y)+N)/windowStep)).*...
    %         buffer(x(N+1:end),N,N-windowStep),NFFT);if dbg, toc;end
    
    % disp(size(repmat(window,1,ceil((length(y)-N)/windowStep))));
    % disp(size(buffer(y,N,N-windowStep,'nodelay')));
    
    [Y ~]=buffer(y,N,N-windowStep,'nodelay');
    Sfft = fft(repmat(window,1,size(Y,2)).*Y,NFFT);
    if dbg, toc;end
    
end



f1 = floor((portion(1)/100)*length(f))+1;
f2 = floor((portion(2)/100)*length(f));
fRange = [f(f1) f(f2)];

% absSpec = 20*log10(abs(Sfft(i1:i2,:))+1);
% angleSpec = angle(Sfft(i1:i2,:));

% if loop
%     if portion(2)<100
%         i1 = floor((portion(1)/100)*size(Sfft,1))+1;
%         i2 = floor((portion(2)/100)*size(Sfft,1));
%         imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))+1));
%     else
%         imagesc(t,fRange,20*log10(abs(Sfft)+1));
%     end
%     xlabel('Time (sec)');ylabel('Frequency (Hz)');
%     set(gca,'YDir','normal');
%     cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');
%     if dbg, toc; end;
% else
%     i1 = floor((portion(1)/200)*size(Sfft,1)+1);
%     i2 = floor((portion(2)/200)*size(Sfft,1));
%     %         imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))));
%     temp=20*log10(abs(Sfft(i1:i2,:))+1);
%     %         temp(temp<1)=1;
%     imagesc(t,fRange,temp);
%     %         imagesc(t,fRange,angle(Sfft(i1:i2,:)));
%     xlabel('Time (sec)');ylabel('Frequency (Hz)');
%     set(gca,'YDir','normal');
%     cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');
%     if dbg, toc; end;
% end





x = zeros(1,length(y)+2*N);
x(N:length(y)+N-1) = y;
start = 11;

Tstep = floor((length(y)+N)/windowStep);
Sautocorr = zeros(N-start+2,Tstep);

count = 1;
for i = 1:windowStep:length(y)+N
    
    Y = xcorr(x(i:i+N-1),N-1,'coeff');
    Y = Y(N:end);
    Sautocorr(:,count) = Y(start-1:end)';
    
    count = count + 1;
    
end
% Sautocorr = Sautocorr*1/max(max(abs(Sautocorr)));

lags = linspace(start,size(Sautocorr,1)-start+1,size(Sautocorr,1))...
    *1000./fs;
screen = get(0,'screensize');
screen(1:2)=screen(3:4)/4;
screen(3:4)=screen(3:4)/2;
figure;set(gcf,'position',screen);
setappdata(gcf, 'SubplotDefaultAxesLocation', [.08 .085 .91 .87]);

subplot(1,2,1)
if loop
    if portion(2)<100
        i1 = floor((portion(1)/100)*size(Sfft,1))+1;
        i2 = floor((portion(2)/100)*size(Sfft,1));
        imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))+1));
    else
        imagesc(t,fRange,20*log10(abs(Sfft)+1));
    end
else
    i1 = floor((portion(1)/200)*size(Sfft,1))+1;
    i2 = floor((portion(2)/200)*size(Sfft,1));
    imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))+1));
end


title('Spectrogram');
xlabel('Time (sec)');ylabel('Frequency (Hz)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');

subplot(1,2,2)
imagesc(t,lags,Sautocorr);colorbar;
title('Autocorrelogram');
xlabel('Time (sec)');ylabel('Lags (ms)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Correlation (-1 to 1)');
drawnow;



drawnow;

