%% cdot
% [dCdt] = cdot(nx, cx)
%
% Integrates one connection

%%
function [dCdt] = cdot(M, tx, nx, cx)

%% Initialize variables and stimulus
% global M
%dbg = 0;

con  = M.n{nx}.con{cx};
C = con.C;
kappa = con.kappa;
lambda = con.lambda;
mu1 = con.mu1;
mu2 = con.mu2;
e = con.e;
nType = con.nType;
no11 = con.no11;

if con.nSourceClass == 1
    y1 = stimulusRun(tx, M.s{con.source});
else
    y1 = M.n{con.source}.z;
    n1 = M.n{con.source};
end
z2   = M.n{con.target}.z;
n2   = M.n{con.target};

%%   Input to connection rule
% $X = \kappa*\vec{z}\vec{z'}^T$

switch nType
    case 1  % 1freq
        X = kappa.*(z2 * y1');
        
    case 2  % 2freq
        Z1 = repmat(y1', n2.N, 1); % conjugate transpose is what we want
        Z2 = repmat(z2 , 1, n1.N);
        N = con.NUM; D = con.DEN;
        X  = kappa.*(e.^((N+D-2)/2)).*((Z2.^D) .* (Z1.^N));
        
    case {3, 4} % 3freq or 3freqall
        Z1 = y1(con.IDX1); Z1(~con.CON1) = conj(Z1(~con.CON1));
        Z2 = y1(con.IDX2); Z2(~con.CON2) = conj(Z2(~con.CON2));
        Z  = z2(con.IDXZ);
        N1 = con.NUM1; N2 = con.NUM2; D = con.DEN;
        Z1N1 = (Z1.^N1);
        Z2N2 = (Z2.^N2);
        ZD  = (Z .^D);
        X = kappa.*(e.^((N1+N2+D-2)/2)).*ZD.*Z1N1.*Z2N2;
        
    case 5  % all2freq
        if no11 % remove 1:1 (and subsequent n:n) monomials
            X = kappa.*(P(e, z2) * P(e, y1') - P(e^2, z2*y1'));
        else
            X = kappa.*(P(e, z2) * P(e, y1'));
        end
        
    case 6  % allfreq
        if no11
            X = kappa.*(P(e, z2) * P_new(e, y1') ...
                - P(e^2, z2*y1') .* repmat(A(e^2, abs(y1').^2), n2.N, 1));
        else
            X = kappa.*(P(e, z2) * P_new(e, y1'));
        end
        
    case 7  % active
        if no11
            X = kappa.*((sqrt(e)*z2.*P(e, z2)) * y1');
        else
            X = kappa.*(P(e, z2) * y1');
        end
end

X = single(X);

%% The differential equation
% $\dot{C} = C (\lambda + \mu_1 C^2+ \epsilon\mu_2 \frac{|C|.^4}{1-\sqrt{\epsilon}*|C|.^2}$ + X

%% oops ...not scaling by frequency
% Why not con.e?
dCdt = C.*(lambda + mu1.*abs(C).^2 + e*mu2.*(abs(C).^4)./(1-e*abs(C).^2)) + X;
end

%% Nonlinear Function Definitions
function y = P(epsilon, x)
y = ( x ./ (1 - sqrt(epsilon)*x) );
end

function y = P_new(epsilon, x)
y = ( x ./ (1 - sqrt(epsilon)*x) ) .* ( 1 ./ (1 - sqrt(epsilon)*conj(x) ));
%y = y - Pc(epsilon, x);
end

function y = A(epsilon, z)
y = ( 1 ./ (1 - sqrt(epsilon)*conj(z) ));
end

function y = Pc(epsilon, x)
y = ( sqrt(epsilon)*x.*conj(x) ./ (1 - epsilon*x.*conj(x)) );
end

function y = Ac(epsilon, x, z)
y = ( sqrt(epsilon)*x.*conj(z) ./ (1 - epsilon*x.*conj(z)) );

end

function y = H(epsilon, r)
y = (epsilon * r.^4 ./ (1- epsilon * r.^2) );
end

function y = Sc(epsilon, z2, z1)
y = ( sqrt(epsilon)*z2*z1' ./ (1 - sqrt(epsilon)*z2*z1') );
% should this be sqrt(epsilon) in the demoninator?
% Ask Ji Chul to verify
end

%% Interpolate stimulus value for given t
function y = stimulusRun(tx, s)

%tx is an index into stimulus data
%check out of bounds. Clamp if needed.

% if tx >= s.lenx %much faster than (length(s.x))
%     y = s.x(ext,s.lenx);
% end

%First check if index is integer or not
rm = tx - floor(tx); %NOTE don't name this 'rem'! That's a func and will slow things!
if rm == 0
    %integer, so valid index
    y = s.x(:,tx);
else
    %We're between indices, so linear interpolation
    y = s.x(:,tx-rm) * (1-rm) + s.x(:,tx+1-rm) * rm;
    %(t-rm) & (t+1-rm) are much faster than floor(t) & ceil(t)
end
end
