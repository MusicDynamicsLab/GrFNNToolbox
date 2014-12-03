function myAutocorr(y,NFFT,fs,varargin)

if nargin<3,error('allMyFreqs needs at least 3 inputs: signal vector, fft size, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
elseif size(y,2)==1
    y=y.';
end
if isempty(varargin)
    portion=[0 100];
    N=2000;
    windowStep=50;
elseif length(varargin)==1
    portion=varargin{1};
    N=2000;
    windowStep=50;
elseif length(varargin)==2
    portion=varargin{1};
    N=varargin{2};
    windowStep=50;
elseif length(varargin)==3
    portion=varargin{1};
    N=varargin{2};
    windowStep=varargin{3};
else
    error('myFreqs only takes 6 inputs')
end
if ~isreal(y)
    warning('Input signal is complex; only real portion taken')
end
y=real(y);

t=0:1/fs:length(y)/fs;

if mod(N,2),N=N+1;end

Tstep = floor((length(y)+N)/windowStep);

x = zeros(1,length(y)+2*N);
x(N:length(y)+N-1) = y;
start = 11;

Sautocorr = zeros(N-start+2,Tstep);

count = 1;
for i = 1:windowStep:length(y)+N
    
    Y = xcorr(x(i:i+N-1),N-1,'coeff');
    Y = Y(N:end);
    Sautocorr(:,count) = Y(start-1:end)';
    
    count = count + 1;
    
end
% Sautocorr = Sautocorr*1/max(max(abs(Sautocorr)));

lags = linspace(start,size(Sautocorr,1),size(Sautocorr,1))...
    *1000./fs;


imagesc(t,lags,Sautocorr);colorbar;
% title('Autocorrelogram');
xlabel('Time (sec)');ylabel('Lags (ms)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Correlation (-1 to 1)');