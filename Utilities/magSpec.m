%% magSpec
%  y = magSpec(x,varargin)
%
%  Computes magnitude spectrum y in dB of input signal x. If there are no
%  other input arguments, the fft length NFFT defaults to the length of x. 
%  If there is a second input argument, it is the fft length NFFT. If there
%  is a third input argument, it is the reference of the dB transform. Else
%  the reference is 1.

function y = magSpec(x,varargin)

if isempty(varargin), NFFT=length(x);else NFFT=varargin{1};end
if ~isreal(x)
    warning('Input signal is complex; only real portion taken')
    x=real(x);
end
if length(varargin)>1
    reference=varargin{2};
else
    reference=1;
end
flag=0;
if size(x,1)==1, x=x';flag=1;end
if length(varargin)>1
    correction=1/size(x,1);
    y=abs(fft(x,NFFT))*2;
else
    correction=1;
    y=abs(fft(x,NFFT));
end
len=size(y,1);
ind=floor(len/2);
if flag
    y=20*log10(y(1:ind)*correction/reference)';
else
    y=20*log10(y(1:ind,:)*correction/reference);
end