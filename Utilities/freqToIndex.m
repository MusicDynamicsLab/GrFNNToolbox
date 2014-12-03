%% freqToIndex
%  FreqToIndex provides the closest index in network n for a given
%  eigenfrequency
%
%  index = FreqToIndex(n, freq)
%
%   x = abs(n.f - freq);
%   index = find(x == min(x));
%

%%
function index = freqToIndex(n, freq)

switch n.fspac
    case 'lin'
        [c,i] = min(abs(n.f - freq));
        index = i;
    case 'log'
        [c,i] = min(abs(log2(n.f / freq)));
        index = i;
    otherwise
        error('Unknown network type')
end

%% revision history
% 1/27/12 JCK and KDL changed method to work for lin and log spacing,
% depending on n.fspac