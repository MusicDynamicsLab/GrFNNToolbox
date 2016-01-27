%% zdot
% [dzdt] = zdot(M, nx, varargin)
%
% Integrates a single network M.n{nx}

%%
function [dzdt] = zdot(M, nx, varargin)
%% Initialize variables and stimulus

% global M

n   = M.n{nx};
z   = n.z;
a   = n.a;
b1  = n.b1;
b2  = n.b2;
e   = n.e;

%% Add input

x = 0;

for cx = 1:length(n.con)
    con = n.con{cx};
    if con.nSourceClass == 1
        y = M.s{con.source}.z;
    else
        y = M.n{con.source}.z;
    end
    
    switch con.nType    % cases ordered by frequency of use
        
        case 1  % 1freq
            x = x + con.w.*(con.C*y);
            
        case 6  % all2freq
            if con.no11
                x = x + con.w .* sum(con.C.*( A(e, z)*P(e, y.') ...
                    - P(e^2, conj(z)*y.')./repmat(conj(z), 1, con.sourceN) ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*P(e, y.') ), 2);
            end
            
        case 7  % allfreq
            if con.no11
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ...
                    - P(e^2, conj(z)*y.').*((1./conj(z))*P(e^2, abs(y.').^2)) ), 2);
            else
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ), 2);
            end
            
        case 2  % 2freq
            NUM = con.NUM;
            DEN = con.DEN;
            Y = repmat(    y.', con.targetN, 1).^ NUM;
            Z = repmat(conj(z) , 1, con.sourceN).^(DEN-1);
            x = x + con.w .* ...
                sum(con.C.*con.epsn.*Y.*Z,2);
            
        case 5  % active
            if con.no11
                x = x + con.w .* sum(con.C.*( (sqrt(e)*conj(z).*A(e, z))*y.' ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*y.' ), 2);
            end
            
        otherwise % 3freq and 3freqall
            Y1 = y(con.IDX1); Y1(con.CON1) = conj(Y1(con.CON1));
            Y2 = y(con.IDX2); Y2(con.CON2) = conj(Y2(con.CON2));
            Z  = conj(z(con.IDXZ));
            NUM1 = con.NUM1;
            NUM2 = con.NUM2;
            DEN = con.DEN;
            x_int = con.C.*con.epsn.*(Y1.^NUM1).*(Y2.^NUM2).*(Z.^(DEN-1));
            % had to conjugate to make it work right
            x = x + con.w .* sum(x_int,2);
            
    end
end

%% The differential equation
% $\dot{\vec{z}} = \vec{z} (a + b_1 \vec{z}^2) +x$
dzdt = z.*(a + b1.*abs(z).^2 + e*b2.*(abs(z).^4)./(1-e*abs(z).^2)) + x;

%%  Nonlinear Function Definitions
function y = P(epsilon, x)
y = ( x ./ (1 - sqrt(epsilon)*x) );

function y = P_new(epsilon, x)
y = ( x ./ (1 - sqrt(epsilon)*x) ) .* ( 1 ./ (1 - sqrt(epsilon)*conj(x) ));
%y = y - Pc(epsilon, x);

function y = A(epsilon, z)
y = ( 1 ./ (1 - sqrt(epsilon)*conj(z) ));

function y = Pc(epsilon, x)
y = ( sqrt(epsilon)*x.*conj(x) ./ (1 - epsilon*x.*conj(x)) );

function y = Ac(epsilon, x, z)
y = ( sqrt(epsilon)*x.*conj(z) ./ (1 - epsilon*x.*conj(z)) );

function y = H(epsilon, r)
y = (epsilon * r.^4 ./ (1- epsilon * r.^2) );
