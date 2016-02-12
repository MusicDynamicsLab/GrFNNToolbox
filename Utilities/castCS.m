%% castCS
%  output = castCS(input)
%
%  Casts input to be output complex and single, regardless of what input 
%  was beforehand. This is required of all numbers in the integrator and 
%  dot functions to make sure we always compute in single precision for 
%  faster speeds. Some C code requires these numbers to be complex by 
%  default as well when interfacing with mex functions. 

function output = castCS(input)

if isreal(input) && isa(input,'double')
    output = complex(single(input));
elseif isreal(input) && isa(input,'single')
    output = complex(input);
elseif ~isreal(input) && isa(input,'double')
    output = single(input);
else
    output = input;
end