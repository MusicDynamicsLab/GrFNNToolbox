function [n1, n2, m] = farey3(x1, x2, z, tol, nvec, nl, mvec, ml, N1, N2, M)

%% Find the sums that equal 0
% See Hoppenstadt & Izhikevich, 1997, pg 173, 265
S = N1*x1 + N2*x2 - M*z;
ix3 = find(abs(S) < tol*z); % & abs(N1)+abs(N2)+abs(M) < 5);

if isempty(ix3)
    n1 = 0; n2 = 0; m = 1;
    return
end

%% Must calculate the row, column & page from ix 
%  because of the way find works (page-major order)
p0 = floor((ix3-1)/(nl*nl)); p = p0+1;
ix2 = ix3 - p0*(nl*nl);
c0 = floor((ix2-1)/nl); c = c0+1;
ix1 = ix2 - c0*nl;
r = ix1;

%% Now find the lowest order monomial and return that
n1 = nvec(r); n2 = nvec(c); m = mvec(p);

E = [n1; n2; m];
ix = find(sum(abs(E)) == min(sum(abs(E))));
ix = ix(1);

n1 = n1(ix); n2 = n2(ix); m = m(ix);
% [n1, n2, m]