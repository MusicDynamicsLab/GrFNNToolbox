function  [X1i, X2i, Zi] = inputShapeOther(num1, num2)

% This function returns the proper indices for making the coupling
% matrices, and their associated exponent matrices
% See test3Freq for a scritped example

%% First, make the basic matrices
x1 = (1:num1)';
x2 = (1:num1);
z  = (1:num2)';

X1 = repmat(x1, 1, num1);
X2 = repmat(x2, num1, 1);

%% Only need upper triangular part (above diagonal)
T1 = triu(X1,1); % T1 is a temporary matrix
T2 = triu(X2,1); % T2 is a temporary matrix

%% Now, for each row, that is output (z), combine the nonzero elements of each
%  column into next row of the matrix;
X1i = [];
X2i = [];
for n = 1:num2

    X1i = [X1i; nonzeros(T1)'];
    X2i = [X2i; nonzeros(T2)'];
    % Note: The number of columns will not be equal to the original number
    % of columns
end

Zi = repmat(z , 1, size(X1i,2));
