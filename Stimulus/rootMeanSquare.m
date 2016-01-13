function RMS = rootMeanSquare(x, dim)
% This function works along the first non-singleton dimension of x
if nargin == 1
    dim = 1;
end
RMS = sqrt(mean(x .* conj(x), dim));