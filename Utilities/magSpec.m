function y=magSpec(x,varargin)

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
    y=abs(fft(x,NFFT))+1;
end
len=size(y,1);
ind=floor(len/2);
if flag
    y=20*log10(y(1:ind)*correction/reference)';
else
    y=20*log10(y(1:ind,:)*correction/reference);
end