%% rootMeanSquare
%  RMS = rootMeanSquare(x, dim)
%
%  Computes root mean square RMS of vector or matrix x
%
%  This function works along the first non-singleton dimension of x
%  if dim is not passed.

function RMS = rootMeanSquare(x, dim)

if nargin == 1
    RMS = sqrt(mean(x .* conj(x)));
else
    RMS = sqrt(mean(x .* conj(x), dim));
end