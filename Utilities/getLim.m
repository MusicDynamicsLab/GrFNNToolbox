%% function: Get axis limit for connection display
% Use it as either
% axLim = getLim(network) or
% axLim = getLim(n.f, n.fspac)

function axLim = getLim(varargin)

if isfield(varargin{1}, 'f')
    n = varargin{1};
    f = n.f;
    fspac = n.fspac;
else
    f = varargin{1};
    fspac = varargin{2};
end
    
if strcmp(fspac, 'log')
    f1 = f(1);
    f2 = f(end);
    N = length(f);
    axMin = ((2*N-1)*f1*(f2/f1)^(-1/(2*(N-1))) + f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
    axMax = (f1*(f2/f1)^(-1/(2*(N-1))) + (2*N-1)*f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
    % Need to do this due to problems in using imagesc with log axis
else
    axMin = f(1);
    axMax = f(end);
end
axLim = [axMin axMax];

end