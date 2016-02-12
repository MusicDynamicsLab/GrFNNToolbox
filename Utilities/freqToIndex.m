%% freqToIndex
%  freqToIndex provides the indices for oscillators in a network, n,
%  whose natural frequencies are closest to given frequencies, freq.
%
%  index = freqToIndex(n, freq)
%

%%
function index = freqToIndex(n, freq)

[NF, FREQ] = meshgrid(n.f, freq);

switch n.nFspac
    case 1 % lin spacing
        [~, index] = min(abs(NF - FREQ), [], 2);
    case 2 % log spacing
        [~, index] = min(abs(log2(NF ./ FREQ)), [], 2);
    otherwise
        error('Unknown network type')
end

%% revision history
% 1/27/12 JCK and KDL changed method to work for lin and log spacing,
% depending on n.fspac