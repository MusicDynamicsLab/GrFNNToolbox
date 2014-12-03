function [N1 N2 M] = inputExponents(f, X1i, X2i, Zi)

Fn1 = f(X1i);
Fn2 = f(X2i);
Fm  = f(Zi );

shape = size(Fn1);

fn1 = Fn1(:);
fn2 = Fn2(:);
fm  = Fm (:);

N1 = zeros(size(fn1));
N2 = zeros(size(fn2));
M  = zeros(size(fm ));

%% Set up the the search space for farey3 as a 3D matrix of all possible monomials
nvec = [(-16:-1), (1:16)]; 
nl = length(nvec);
mvec = [ 1:5]; 
ml = length(mvec);

N1in = repmat(nvec', [ 1, nl, ml]);
N2in = repmat(nvec , [nl,  1, ml]);
Min = cat(3, [1*ones(nl,nl), 2*ones(nl,nl), 3*ones(nl,nl), 4*ones(nl,nl), 5*ones(nl,nl)]);
Min = reshape(Min, nl, nl, ml);

%% Not sure how to build 3rd dimension (M) properly
%  ... but this works well enough for now

%% This loop takes lots of time!
N1 = zeros(size(fm));
N2 = zeros(size(fm));
M  = zeros(size(fm));
for fi = 1:length(fm)
    [n1 n2 m] = farey3(fn1(fi), fn2(fi), fm(fi), .02, nvec, nl, mvec, ml, N1in, N2in, Min);
    N1(fi) = n1;
    N2(fi) = n2;
    M (fi) = m ;
    
end

N1 = reshape(N1, shape);
N2 = reshape(N2, shape);
M  = reshape(M , shape);
       