%% zdot
%  dzdt = zdot(M, nx)
%
%  Calculates time derivative for a single network M.n{nx}
%
%  Input arguments:
%  M        Model
%  nx       Network id
%
%  Output:
%  dzdt     Time derivative of oscillator states (column vector)
%

%%
function [dzdt] = zdot(M, nx)
%% Initialize variables and stimulus

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
                    - A(e^2, z*y').*repmat(y.', con.targetN, 1) ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*P(e, y.') ), 2);
            end
            
        case 7  % allfreq
            if con.no11
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ...
                    - A(e^2, z*y').*repmat(y.'.*A(e^2, abs(y.').^2), con.targetN, 1) ), 2);
            else
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ), 2);
            end
            
        case 2  % 2freq
            NUM = con.NUM;
            DEN = con.DEN;
            Y = repmat(sqrt(e)*y.', con.targetN, 1).^ NUM;
            Z = repmat(sqrt(e)*conj(z) , 1, con.sourceN).^(DEN-1);
            x = x + con.w .* sum(con.C.*Y.*Z,2)/sqrt(e);
            
        case 5  % active
            if con.no11
                x = x + con.w .* sum(con.C.*( (sqrt(e)*conj(z).*A(e, z))*y.' ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*y.' ), 2);
            end
            
        otherwise % 3freq and 3freqall
            Y1 = sqrt(e)*y(con.IDX1); Y1(con.CON1) = conj(Y1(con.CON1));
            Y2 = sqrt(e)*y(con.IDX2); Y2(con.CON2) = conj(Y2(con.CON2));
            Z  = sqrt(e)*conj(z(con.IDXZ));
            NUM1 = con.NUM1;
            NUM2 = con.NUM2;
            DEN = con.DEN;
            x_int = con.C.*(Y1.^NUM1).*(Y2.^NUM2).*(Z.^(DEN-1))/sqrt(e);
            x = x + con.w .* sum(x_int,2);
            
    end
end

%% The differential equation
% $\dot{z} = z \left( \alpha + \textrm{i}\omega + (\beta_1 + \textrm{i}\delta_1) |z|^2 + \frac{\epsilon (\beta_2 + \textrm{i}\delta_2) |z|^4}{1-\epsilon |z|^2} \right) + x$
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
