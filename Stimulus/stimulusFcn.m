%% STIMULUSFCN
%  Generates time series for stimulus signal.
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
        linear = 0;
        if length(temp{i}) > 1
            if length(temp{i}) == 2, linear = 1;end
            temp{i} = interp1(linspace(0,t(end),length(temp{i})),temp{i},t);
            if i < 4
                deriv = diff(temp{i})*fs;
                if linear
                    deriv = -[deriv deriv(end)];
                else
                    deriv = -interp1(linspace(deriv(1),deriv(end),...
                        length(deriv)),deriv,linspace(deriv(1),deriv(end),...
                        length(deriv)+1));
                end
                temp{i} = (temp{i} + t.*deriv/2);
            end
        end
    end
    
    fc = temp{1};fm = temp{2};fAM = temp{3};
    ac = temp{4};am = temp{5};aAM = temp{6};
    
%% Create frequency modulation
    
    if isfield(s, 'FM')
        
        if fm ~= 0
            tflr = floor(t.*fm)./fm;
            trem = rem(t,1./fm);
            q = 2*pi*fm.*t;
            
            switch lower(s.FM{a,b}(1:3))
                case 'cos'
                    mi= ((q-sin(q))/2)/2/pi./fm; % integral of -cos(q)
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
    
    if isfield(s, 'AM');
        
        Cfreq = s.CfreqsAM(a,b);
        AM = zeros(size(t));
        
        switch lower(s.AM{a,b})
            case 'cos'
                AM = (Cfreq + aAM.*cos(2*pi*fAM.*t));
            case 'sin'
                AM = (Cfreq + aAM.*sin(2*pi*fAM.*t));
            case 'squ'
                AM = (Cfreq + aAM.*square(2*pi*fAM.*t));
            case 'saw'
                AM = (Cfreq + aAM.*sawtooth(2*pi*fAM.*t));
            otherwise
                AM = 1;
                warning('Unknown AM modulator type, not modulating');
        end
    else
        AM = 1;
    end
    
%% Create carrier and apply modulators
    
    xb = zeros(size(t));
    
    switch lower(s.carrier{a,b}(1:3))
        case 'cos'
            xb = cos(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b)).*AM;
        case 'sin'
            xb = sin(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b)).*AM;
        case 'exp'
            xb = exp(1i*(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b))).*AM;
        case 'saw'
            xb = sawtooth(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b)).*AM;
        case 'squ'
            xb = square(2*pi*fc.*t + 2*pi*fc.*am.*mi + s.Th(a,b)).*AM;
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
