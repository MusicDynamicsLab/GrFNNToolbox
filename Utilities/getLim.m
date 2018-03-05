%% getLim
%  axLim = getLim(varargin)
%
%  Gets axis limit to be used for imagesc of connection. Needed to solve
%  problems in using imagesc with log axis.
%
%  Use it as either
%   axLim = getLim(network) or
%   axLim = getLim(network.f, network.nFspac)
%     nFspac is 1 for linear spacing or 2 for log
%

%%
function axLim = getLim(varargin)

if isfield(varargin{1}, 'f')
    n = varargin{1};
    f = n.f;
    nFspac = n.nFspac;
else
    f = varargin{1};
    nFspac = varargin{2};
end
    
switch nFspac
    case 2 % log spacing
        f1 = f(1);
        f2 = f(end);
        N = length(f);
        axMin = ((2*N-1)*f1*(f2/f1)^(-1/(2*(N-1))) + f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
        axMax = (f1*(f2/f1)^(-1/(2*(N-1))) + (2*N-1)*f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
        % Need to do this due to problems in using imagesc with log axis
    case 1 % lin spacing
        axMin = f(1);
        axMax = f(end);
end
axLim = [axMin axMax];

end