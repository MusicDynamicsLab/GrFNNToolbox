%% stimulusFcn
%  x = stimulusFcn(t, s, a)
%  Generates analytical time series x for stimulus signal given t, 
%  stimulus time vector, s, stimulus structure, and a, section of 
%  stimulus (row number of s.ts)
%
%  Called by stimulusMake; not a standalone function.

%%
function x = stimulusFcn(t, s, a)


x = zeros(size(t));               % Initialize signal vector
fs = s.fs;

%% Loop through number of components of signal

for b = 1:length(s.fc(a,:))       % For each component of the section
    
    fc = s.fc{a,b};ac = s.ac{a,b};
    
    if isfield(s, 'FM'), fm = s.fFM{a,b};am = s.aFM{a,b};
    else
        fm = 0;am = 0;
    end
    
    if isfield(s, 'AM'), fAM = s.fAM{a,b};aAM = s.aAM{a,b};
    else
        fAM = 0;aAM = 0;
    end
    
%% Create complete amplitude/frequency vectors if time-varying
    
    temp = {fc fm fAM ac am aAM};
    
    for i = 1:6
        if length(temp{i}) > 1
            temp{i} = interp1(linspace(0,t(end),length(temp{i})),temp{i},t);
        else
	        temp{i} = repmat(temp{i},1,length(t));
        end
    end
    
    fc  = temp{1};
    fm  = temp{2};
    fAM = temp{3};
    ac  = temp{4};
    am  = temp{5};
    aAM = temp{6};
    
%% Create frequency modulation
    
    if isfield(s, 'FM')
        
        if any(fm ~= 0)
            tflr = floor(t.*fm)./fm;
            trem = rem(t,1./fm);
%             q = 2*pi*fm.*t;
            q = 2*pi*(cumsum(fm)-fm(1))/fs;
            
            switch lower(s.FM{a,b}(1:3))
                case 'sin'
%                     mi= ((q-sin(q))/2)/2/pi./fm; % integral of -cos(q)
                    moddedFreqs = am.*sin(q);
                    moddedFreqs = cumsum(moddedFreqs)-moddedFreqs(1);
                    mi= moddedFreqs/fs;
                case 'cos'
                    moddedFreqs = am.*cos(q);
                    moddedFreqs = cumsum(moddedFreqs)-moddedFreqs(1);
                    mi= moddedFreqs/fs;
                case 'saw'
                    mi= tflr/2 + (fm.*trem.^2)/2; % (fm*t.^2)/2 is integral of fm*t
                case 'squ'
                    hlin = zeros(size(trem)); ii = find(trem>1./fm/2); % yuk!!!
                    hlin(ii) = trem(ii)-1./fm/2;
                    mi= (tflr + 2*hlin)/2; % (2*t)/2 is integral of 1
                otherwise
                    am = 0;
                    mi = 0;
                    warning('Unknown FM modulator type, not modulating');
            end
        else
            am = 0;
            mi = 0;
        end
    else
        am = 0;
        mi = 0;
    end
    
%% Create amplitude modulation
    
    if isfield(s, 'AM')
        
        Cfreq = s.CfreqsAM(a,b);
        q = 2*pi*(cumsum(fAM)-fAM(1))/fs;
        
        switch lower(s.AM{a,b})
            case 'cos'
                AM = (Cfreq + aAM.*cos(q));
            case 'sin'
                AM = (Cfreq + aAM.*sin(q));
            case 'squ'
                AM = (Cfreq + aAM.*square(q));
            case 'saw'
                AM = (Cfreq + aAM.*sawtooth(q));
            otherwise
                AM = 1;
                warning('Unknown AM modulator type, not modulating');
        end
    else
        AM = 1;
    end
    
%% Create carrier and apply modulators
    
    xb = zeros(size(t));
    freqsIntegral = 2*pi*(cumsum(fc)-fc(1))/fs;
    
    switch lower(s.carrier{a,b}(1:3))
        case 'cos'
            xb = cos(freqsIntegral + 2*pi*fc.*mi + s.Th(a,b)).*AM;
        case 'sin'
            xb = sin(freqsIntegral + 2*pi*fc.*mi + s.Th(a,b)).*AM;
        case 'exp'
            xb = exp(1i*(freqsIntegral + 2*pi*fc.*mi + s.Th(a,b))).*AM;
        case 'saw'
            xb = sawtooth(freqsIntegral + 2*pi*fc.*mi + s.Th(a,b)).*AM;
        case 'squ'
            xb = square(freqsIntegral + 2*pi*fc.*mi + s.Th(a,b)).*AM;
        case 'bmp'
            x1 = sqrt(.75)*cos(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b));
            xb = (x1./(1-x1)-1).*AM;
        case 'noi'
            for i = 1:length(t)
                xb(i) = 2*rand-1;
            end
            xb = xb.*AM;
        case 'pls' %DH added for use with midi
            xb = exp(4*sin(2*pi*fc*t))/exp(4.05);
            s.sc = -1;
        otherwise
            error('Unknown carrier wave type');
    end
    
    y = ac.*xb;
    
%% Apply delay-and-add circuit
    
    if isfield(s, 'iter') && s.iter(a,b) ~= 0 && s.Niter(a,b) ~= 0
        iter = floor(s.iter(a,b) * fs);
        Niter = s.Niter(a,b);
        y = stimulusIter(y, iter, Niter);
    end
    
%% Add component to signal
    
    x = x + y;
    
end
