%% cdot
% [dCdt] = cdot(nx, cx)
%
% Integrates one connection

%%
function [dCdt] = cdot(nx, cx)

%% Initialize variables and stimulus
global M
%dbg = 0;

con  = M.n{nx}.con{cx};
C = con.C;
kappa = con.kappa;
lambda = con.lambda;
mu1 = con.mu1;
mu2 = con.mu2;
e = con.e;
type = con.type;
no11 = con.no11;

n1   = M.n{con.n1};
n2   = M.n{con.n2};

%%   Input to connection rule
% $X = \kappa*\vec{z}\vec{z'}^T$

if strcmpi(type, '1freq')
    X = kappa.*(n2.z * n1.z');
elseif strcmpi(type, '2freq')
    Z1 = repmat(n1.z', n2.N, 1); % conjugate transpose is what we want
    Z2 = repmat(n2.z , 1, n1.N);
    N = con.N; D = con.D;
    X  = kappa.*(e.^((N+D-2)/2)).*((Z2.^D) .* (Z1.^N));
elseif strcmpi(type, '3freq')
    Z1 = n1.z(con.IDX1); Z1(~con.CON1) = conj(Z1(~con.CON1));
    Z2 = n1.z(con.IDX2); Z2(~con.CON2) = conj(Z2(~con.CON2));
    Z  = n2.z(con.IDZ);
    N1 = con.NUM1; N2 = con.NUM2; D2 = con.DEN2;
    
    Z1N1 = (Z1.^N1);
    Z2N2 = (Z2.^N2);
    ZD2  = (Z .^D2);

    X = kappa.*(e.^((N1+N2+D2-2)/2)).*ZD2.*Z1N1.*Z2N2;
elseif strcmpi(type, 'All2freq')
    if no11 % remove 1:1 (and subsequent n:n) monomials
        X = kappa.*(P(e, n2.z) * P(e, n1.z') - P(e^2, n2.z*n1.z'));
    else
        X = kappa.*(P(e, n2.z) * P(e, n1.z'));
    end
elseif strcmpi(type, 'Allfreq')
    if no11
        X = kappa.*(P(e, n2.z) * P_new(e, n1.z') ...
            - P(e^2, n2.z*n1.z') .* repmat(A(e^2, abs(n1.z').^2), n2.N, 1));
    else
        X = kappa.*(P(e, n2.z) * P_new(e, n1.z'));
    end
end

X = single(X);

%% The differential equation
% $\dot{C} = C (\lambda + \mu_1 C^2+ \epsilon\mu_2 \frac{|C|.^4}{1-\sqrt{\epsilon}*|C|.^2}$ + X

%% oops ...not scaling by frequency
% Why not con.e?
dCdt = C.*(lambda + mu1.*abs(C).^2 + e*mu2.*(abs(C).^4)./(1-sqrt(e)*abs(C).^2)) + X;
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
