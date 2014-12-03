%% IndexToFreq
% FreqToIndex provides the eigenfrequency of an oscillator in a network n
% given its index
% Note: can also just query array 'n.f' directly

%%
function freq = IndexToFreq(n, idx)


freq = n.f(idx);
