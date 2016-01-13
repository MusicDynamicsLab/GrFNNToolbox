%% zdot
% [dzdt] = zdot(t, nx, varargin)
%
% Integrates a single network

%%
function [dzdt] = zdot(M, tx, nx, varargin)
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
        y1 = stimulusRun(tx, M.s{con.source});
    else
        y1 = M.n{con.source}.z;
    end
    
    switch con.nType
        case 1  % 1freq
            x = x + con.w.*(con.C*y1);
            
        case 2  % 2freq
            N1 = M.n{con.source}.N;
            N2 = n.N;
            N = con.NUM;
            D = con.DEN;
            Y1 = repmat(    y1.', N2, 1).^ N;
            Z2 = repmat(conj(z) , 1, N1).^(D-1);
            x = x + con.w .* ...
                sum(con.C.*(e.^((N+D-2)/2)).*Y1.*Z2,2);
            
        case {3, 4} % 3freq and 3freqall
            Y1 = y1(con.IDX1); Y1(con.CON1) = conj(Y1(con.CON1));
            Y2 = y1(con.IDX2); Y2(con.CON2) = conj(Y2(con.CON2));
            Z  = conj(z(con.IDXZ));
            N1 = con.NUM1;
            N2 = con.NUM2;
            D = con.DEN;
            x_int = con.C.*(e.^((N1+N2+D-2)/2)) ...
                .*(Y1.^N1).*(Y2.^N2).*(Z.^(D-1));
            % had to conjugate to make it work right
            x = x + con.w .* sum(x_int,2);
            
        case 5  % all2freq
            if con.no11
                if con.nSourceClass == 1
                    N1 = M.s{con.source}.Nch;
                else
                    N1 = M.n{con.source}.N;
                end
                x = x + con.w .* sum(con.C.*( A(e, z)*P(e, y1.') ...
                    - P(e^2, conj(z)*y1.')./repmat(conj(z),1,N1) ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*P(e, y1.') ), 2);
            end
            
        case 6  % allfreq
            if con.no11
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y1.') ...
                    - P(e^2, conj(z)*y1.').*((1./conj(z))*P(e^2, abs(y1.').^2)) ), 2);
            else
                x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y1.') ), 2);
            end
            
        case 7  % active
            if con.no11
                x = x + con.w .* sum(con.C.*( (sqrt(e)*conj(z).*A(e, z))*y1.' ), 2);
            else
                x = x + con.w .* sum(con.C.*( A(e, z)*y1.' ), 2);
            end
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
