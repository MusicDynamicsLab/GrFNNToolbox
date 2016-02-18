%% mdlAutocorr
%  mdlAutocorr(y,fs,varargin)
%
%  Plotting function. Produces an autocorrelogram of signal vector y. Input
%  argument for sampling frequency fs is required. Optional arguments must
%  come in order afterwards, including N and windowStep. N is the desired
%  window length in samples, and windowStep is the desired window step size
%  in samples.

%%
function mdlAutocorr(y,fs,varargin)

if nargin<2,error('mdlAutocorr needs at least 2 inputs: signal vector, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
end

y=y(:);

if isempty(varargin)
    N=ceil(length(y)^(1/1.4));
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==1
    N=varargin{1};
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==2
    N=varargin{1};
    windowStep=varargin{2};
else
    error('mdlAutocorr only takes 4 inputs')
end
if ~isreal(y)
    warning('Input signal is complex; only real portion taken')
    y=real(y);
end

t=(0:length(y))/fs;

% if mod(N,2),N=N+1;end

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

lags = linspace(start,size(Sautocorr,1)-start+1,size(Sautocorr,1))...
    *1000./fs;

imagesc(t,lags,Sautocorr);colorbar;
title('Autocorrelogram');
xlabel('Time (sec)');ylabel('Lags (ms)');
set(gca,'YDir','normal');
cbar = colorbar;set(get(cbar,'ylabel'),'string','Correlation (-1 to 1)');