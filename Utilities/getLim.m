%% function: Get axis limit for connection display
function axLim = getLim(n)
if strcmp(n.fspac, 'log')
    f1 = n.f(1);
    f2 = n.f(end);
    N = n.N;
    axMin = ((2*N-1)*f1*(f2/f1)^(-1/(2*(N-1))) + f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
    axMax = (f1*(f2/f1)^(-1/(2*(N-1))) + (2*N-1)*f2*(f2/f1)^(1/(2*(N-1))))/(2*N);
    % Need to do this due to problems in using imagesc with log axis
else
    axMin = n.f(1);
    axMax = n.f(end);
end
axLim = [axMin axMax];

end