%% fareyratio
%  [num, denom, level, error] = fareyratio(fractions, tol)
%
%  Obtains simplest integer ratios, in the form of numerators and
%  denominators, for given fractions. It searches through the Farey
%  sequences of increasing order until an integer ratio is found within
%  tolerance range. See https://en.wikipedia.org/wiki/Farey_sequence
%
%  Input arguments:
%  fractions        A vector of fractions
%  tol              Tolerance to be met to choose an integer ratio
%                   near a given fraction (default is .01, i.e. 1%)
%
%  Ouput:
%  num              Numerators for integer ratios
%  denom            Denominators for integer ratios
%  level            Orders of integer ratios in Farey tree (low level means
%                   small numerator and denominator)
%  error            Difference between given fractions and obtained integer
%                   ratios
%

%%
function [num, denom, level, error] = fareyratio(fractions, tol)

if nargin < 2
    tol = .01;
end

num = [];
denom = [];
level = [];
error = [];

for f = fractions
    
    if f > 1
        frac = 1/f;
    else
        frac = f;
    end
    
    ln = 0;
    ld = 1;
    
    rn = 1;
    rd = 1;
    
    l = 1;
    
    if (abs(frac - ln/ld) <= frac*tol)
        n = ln;
        d = ld;
        e = abs(frac - ln/ld);
        
    elseif (abs(frac - rn/rd) <= frac*tol)
        n = rn;
        d = rd;
        e = abs(frac - rn/rd);
        
    else
        cn = ln+rn;
        cd = ld+rd;
        l  = l + 1;
        
        while (abs(frac - cn/cd) > frac*tol)
            
            if frac > cn/cd
                ln=cn;ld=cd;
            else
                rn=cn;rd=cd;
            end
            
            cn = ln+rn;
            cd = ld+rd;
            l  = l + 1;
            
        end
        
        n = cn;
        d = cd;
        e = abs(frac - cn/cd);
    end
    
    if f > 1
        tmp = n;
        n = d;
        d = tmp;
    end
    %   e = abs(f - n/d);
    
    num   = [num   n];
    denom = [denom d];
    level = [level l];
    error = [error e];
end
