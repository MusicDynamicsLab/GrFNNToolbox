function RMS = rms(x)
% This function works along the first non-singleton dimension of x
RMS = sqrt(mean(x .* conj(x)));