%% STIMULUSRAMP
%  Ramps signal up from zero and down back to zero, acc. to linear or nonlinear scale
%  If y is a matrix, stimulusRamp works along columns.
%  So,
%  x = stimulusRamp(signal, time in sec of ramps, exponent of scale (+), sampling rate of signal)
%
%  Examples:
%  S = stimulusRamp(x, 0.05, 1, 44100);
%  Ramps the start up and the end down of signal x over the course of 50 ms each, linearly.
%  S = stimulusRamp(x, 0.05, 4, 44100);
%  Does the same thing except ramps are more sudden (by factor of 4) than linear. A number between
%  0 and 1 here makes ramps less sudden. No negative values.

%%
function x = stimulusRamp(y, c, p, fs)


if nargin ~= 4
    error('stimulusRamp takes four arguments')
elseif isempty(y)
    error('Signal is empty')
end

cutoff = round(fs*c);

row = 0;
maxSamp = cutoff;
if size(y,1) == 1 && size(y,2) > 1
    row = 1;
    y = y.';
end

len = size(y,1);

scfcn = ones(len,1);

if 2 * cutoff == len
    warning('Stimulus section length equal to summed ramp times')
elseif 2 * cutoff > len && len > cutoff
    warning('Stimulus section length shorter than summed ramp times')
    maxSamp = len / 2;
elseif len == cutoff
    error('Stimulus section length equal to ramp time')
elseif len < cutoff
    error('Stimulus section length shorter than ramp time')
end


for i = 1:maxSamp
    scfcn(i) = ((i-1)/cutoff)^p;
    scfcn(len+1-i) = ((i-1)/cutoff)^p;
end

scfcn = repmat(scfcn,1,size(y,2));

if row
    x = y .* scfcn;
    x = x.';
else
    x = y .* scfcn;
end
