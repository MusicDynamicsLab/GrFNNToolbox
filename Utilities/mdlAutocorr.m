%% mdlAutocorr
%  mdlAutocorr(y,fs,varargin)
%
%  Plotting function. Produces an autocorrelogram of signal vector y. Input
%  argument for sampling frequency fs is required. Optional arguments must
%  come in order afterwards, including portion, windowSize, and windowStep.
%  portion is a two-element vector of percentages indicating a window of
%  lag amounts, from zero to the maximum allowable given length(y), to be
%  displayed. windowSize is the desired window length in samples, and
%  windowStep is the desired window step size in samples.

%%
function [h,cbar] = mdlAutocorr(y,fs,varargin)

if nargin<2,error('mdlAutocorr needs at least 2 inputs: signal vector, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
end

y=y(:);

if isempty(varargin)
    portion=[0 100];
    windowSize=ceil(length(y)^(1/1.4));
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==1
    portion=varargin{1};
    windowSize=ceil(length(y)^(1/1.4));
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==2
    portion=varargin{1};
    windowSize=varargin{2};
    windowStep=ceil(length(y)^(1/2.5));
elseif length(varargin)==3
    portion=varargin{1};
    windowSize=varargin{2};
    windowStep=varargin{3};
else
    error('mdlAutocorr only takes 5 inputs')
end
if ~isreal(y)
    warning('Input signal is complex; only real portion taken')
    y=real(y);
end

t=(0:length(y))/fs;

Tstep = floor((length(y)+windowSize)/windowStep);  % Calculate the size of the autocorrelogram matrix in the time dimension

x = zeros(1,length(y)+2*windowSize);  % Create vector of zeros size of input plus a windowSize on both sides
x(windowSize:length(y)+windowSize-1) = y;  % Write input into the middle of that

Sautocorr = zeros(windowSize,Tstep);  % Autocorrelogram matrix will have size in lags dimension same as input

count = 1;
for i = 1:windowStep:length(y)+windowSize
    
    Y = xcorr(x(i:i+windowSize-1),'coeff');  % Calculate autocorrelation for each window of the zero-padded vector
    Y = Y(windowSize:end);  % and only take the last half

    Sautocorr(:,count) = Y';
    
    count = count + 1;
    
end

lags = (0:windowSize)*1000/fs;

inds=round(portion/100*size(Sautocorr,1));

if inds(1)==0, inds(1)=1;end  % stupid but simple

h = imagesc(t,lags(inds(1):inds(2)),Sautocorr(inds(1):inds(2),:));
title('Autocorrelogram')
xlabel('Time (sec)')
ylabel('Lags (ms)')
set(gca,'YDir','normal')
cbar = colorbar;
set(get(cbar,'ylabel'),'string','Correlation (-1 to 1)')
