function  [X1i, X2i, Zi] = inputShapeInternal(N)

% This function returns the proper indices for making the coupling
% matrices, and their associated exponent matrices
% See test3Freq for a scritped example

%% First, make the basic matrices
x1 = (1:N)';
x2 = (1:N);
z  = (1:N)';

X1 = repmat(x1, 1, N);
X2 = repmat(x2, N, 1);

%% Only need upper triangular part (above diagonal)
T1 = triu(X1,1); % T1 is a temporary matrix
T2 = triu(X2,1); % T2 is a temporary matrix

%% Now, for each row, eliminate row & col corresponding to that row number
%  This eliminates self-connections, the resulting rows will correspond to the
%  subharmonic (the Zi)
% This only really works for internal connections, because it uses index as
% a proxy for frequency
% Also, It may not even be necessary, cause in the next step, these will
% get zero exponents
X1i = [];
X2i = [];
for n = 1:N
    ix = setdiff((1:N), n);
    T1p = T1(ix,ix);
    T2p = T2(ix,ix);
    % Combine the nonzero elements of each column into next row of the
    % matrix; column 1 has no 1's in it, column 2 has no 2's in it, etc
    X1i = [X1i; nonzeros(T1p)'];
    X2i = [X2i; nonzeros(T2p)'];
    % Note: The number of columns will not be equal to the original number
    % of columns
end

Zi = repmat(z , 1, size(X1i,2));
