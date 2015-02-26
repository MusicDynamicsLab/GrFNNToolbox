function struct = myFreqs(y,NFFT,fs,varargin)

if nargin<3,error('myFreqs needs at least 3 inputs: signal vector, fft size, and samp. freq.');end
if size(y,1)>1 && size(y,2)>1
    error('Signal must be only a vector');
end
y=y(:);
if nargout>1, error('Can have only one arg out: Struct with .spec (spectrogram mat), .t (time vector) and .f (freq vector)');end
loop=0;
dbg=1;
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
%     Sfft = fft(repmat(window,1,size(Y,2)).*Y,NFFT);
    Sfft = fft(bsxfun(@times,Y,window),NFFT);
    if dbg, toc;end
    
end



f1 = floor((portion(1)/100)*length(f))+1;
f2 = floor((portion(2)/100)*length(f));
fRange = [f(f1) f(f2)];

% absSpec = 20*log10(abs(Sfft(i1:i2,:))+1);
% angleSpec = angle(Sfft(i1:i2,:));
if nargout==0
    if loop
        if portion(2)<100
            i1 = floor((portion(1)/100)*size(Sfft,1))+1;
            i2 = floor((portion(2)/100)*size(Sfft,1));
            imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))+1));
        else
            imagesc(t,fRange,20*log10(abs(Sfft)+1));
        end
        xlabel('Time (sec)');ylabel('Frequency (Hz)');
        set(gca,'YDir','normal');
        cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');
        if dbg, toc; end;
    else
        i1 = floor((portion(1)/200)*size(Sfft,1)+1);
        i2 = floor((portion(2)/200)*size(Sfft,1));
        %         imagesc(t,fRange,20*log10(abs(Sfft(i1:i2,:))));
        temp=20*log10(abs(Sfft(i1:i2,:))+1);
%         temp(temp<1)=1;
        imagesc(t,fRange,temp);
        %         imagesc(t,fRange,angle(Sfft(i1:i2,:)));
        xlabel('Time (sec)');ylabel('Frequency (Hz)');
        set(gca,'YDir','normal');
        cbar = colorbar;set(get(cbar,'ylabel'),'string','Amplitude (dB)');
        if dbg, toc; end;
    end
else
    if loop
        if portion(2)<100
            i1 = floor((portion(1)/100)*size(Sfft,1))+1;
            i2 = floor((portion(2)/100)*size(Sfft,1));
            struct.spec = 20*log10(abs(Sfft(i1:i2,:))+1);
            struct.t = t;
            struct.f = linspace(fRange(1),fRange(2),size(struct.spec,1));
        else
            struct.spec = 20*log10(abs(Sfft)+1);
            struct.t = t;
            struct.f = linspace(fRange(1),fRange(2),size(struct.spec,1));
        end
        if dbg, toc; end;
    else
        i1 = floor((portion(1)/200)*size(Sfft,1))+1;
        i2 = floor((portion(2)/200)*size(Sfft,1));
        struct.spec = 20*log10(abs(Sfft(i1:i2,:))+1);
        struct.t = t;
        struct.f = linspace(fRange(1),fRange(2),size(struct.spec,1));
        if dbg, toc; end;
    end
end
