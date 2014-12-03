%% moduRhythm
%
%   [Nm, Np] = moduRhythm(N, bper, mdep, mfun, mper)
%
%  N = notematrix to be modulated
%  bper = the beat period (1/tempo in Hz)
%  mdep = modulation depth
%  mfun = modulation funtion
%  mper = modulation period in beats

%%
function [Nm, Np, m, CH] = moduRhythm(N, bper, mdep, mfun, mper)


if nargin < 5; mper = 16;   end;
if nargin < 4; mfun = 'sin'; end;
if nargin < 3; mdep =  .2;   end;
if nargin < 2; bper =  .5;   end;
    
b0   = N(1,1);
b    = N(:,1);
db   = N(:,2);
eb1  = [b; b+db];       % events, both onsets and offsets
[eb ix]  = sort(eb1);   % put in order of occurence, remember how order was
                        % changed so we can make them durations again later

ieib = diff(eb);        % inter event interval in beats 

m   = bper*mdep * modfun((eb-b0)/mper, mfun);
m1  = m; m1(2:end)=m1(2:end).*ieib; % sort of funny, but works right
et  = eb .* bper + cumsum(m1);

ieit = diff(et);
[tmp ixr] = sort(ix);   % put back in original unsorted order 
et1 = et(ixr);          % (onset times, offset times)

t  = et1(1:end/2);
dt = et1(end/2+1:end)-t;

CH = [eb1'; eb'; ix'; ixr'; et'; et1'];

Nm = N;
Nm(:,6) = t-min(t);
Nm(:,7) = dt;

Np = Nm;
Np(:,1) = Nm(:,6); % 120 bpm = 2 Hz
Np(:,2) = Nm(:,7);

%% modfun
function m = modfun(t, type)

switch type
    case 'linear'
        m = t;
    case 'sin'
        m = sin(2*pi*t);
    case 'cos'
        m = cos(2*pi*t);
    case 'sawtooth'
        m = sawtooth(2*pi*t);
    case 'square'
        m = square(2*pi*t);
end