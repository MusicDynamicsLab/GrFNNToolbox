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