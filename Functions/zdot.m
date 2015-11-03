%% zdot
% [dzdt] = zdot(t, nx, varargin)
% 
% Integrates a single network

%%
function [dzdt] = zdot(M, stimulus, t, nx, varargin)
%% Initialize variables and stimulus

% global M

n   = M.n{nx};
z   = n.z;
a   = n.a;
b1  = n.b1;
b2  = n.b2;
e   = n.e;
ext = n.ext;
con = n.con;

%%   External stimulus
x = 0;
if ext
  x = stimulusRun(t, stimulus, ext);  % External signal, scalar
end
switch lower(stimulus.inputType)
    case 'allfreq'
        x = n.w .* P_new(e, x) .* A(e, z);
    case 'all2freq'
        x = n.w .* P(e, x) .* A(e, z);
    case 'active'
        x = n.w .* x .* A(e, z);
    otherwise
        x = n.w .* x;
end

%% Connection input
for cx = 1:length(con)
    z1 = M.n{con{cx}.n1}.z;
    if strcmpi(con{cx}.type, '1freq')
        x = x + con{cx}.w.*(con{cx}.C*z1);
    elseif strcmpi(con{cx}.type, '2freq')
        N1 = M.n{con{cx}.n1}.N;
        N2 = n.N;
        N = con{cx}.NUM;
        D = con{cx}.DEN;
        Z1 = repmat(    z1.', N2, 1).^ N;
        Z2 = repmat(conj(z) , 1, N1).^(D-1);
        x = x + con{cx}.w .* ...
            sum(con{cx}.C.*(e.^((N+D-2)/2)).*Z1.*Z2,2); 
                                 % sum along dimension 2
    elseif strcmpi(con{cx}.type(1:5), '3freq')
        Z1 = z1(con{cx}.IDX1); Z1(con{cx}.CON1) = conj(Z1(con{cx}.CON1));
        Z2 = z1(con{cx}.IDX2); Z2(con{cx}.CON2) = conj(Z2(con{cx}.CON2));
        Z  = conj(z(con{cx}.IDXZ));
        N1 = con{cx}.NUM1; 
        N2 = con{cx}.NUM2; 
        D = con{cx}.DEN;

        x_int = con{cx}.C.*(e.^((N1+N2+D-2)/2)) ...
                         .*(Z1.^N1).*(Z2.^N2).*(Z.^(D-1));
        % had to conjugate to make it work right
        x = x + con{cx}.w .* sum(x_int,2);
    elseif strcmpi(con{cx}.type, 'All2freq')
        if con{cx}.no11
            N1 = M.n{con{cx}.n1}.N;
            x = x + con{cx}.w .* sum(con{cx}.C.*( A(e, z)*P(e, z1.') ...
                - P(e^2, conj(z)*z1.')./repmat(conj(z),1,N1) ), 2);
        else
            x = x + con{cx}.w .* sum(con{cx}.C.*( A(e, z)*P(e, z1.') ), 2);
        end
    elseif strcmpi(con{cx}.type, 'Allfreq')
        if con{cx}.no11
            x =  x + con{cx}.w .* sum(con{cx}.C.*( A(e, z)*P_new(e, z1.') ...
                - P(e^2, conj(z)*z1.').*((1./conj(z))*P(e^2, abs(z1.').^2)) ), 2);
        else
            x =  x + con{cx}.w .* sum(con{cx}.C.*( A(e, z)*P_new(e, z1.') ), 2);
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
function y = stimulusRun(t, s, ext)

%t is an index into stimulus data
%check out of bounds. Clamp if needed.
  
% if t >= s.lenx %much faster than (length(s.x))
%     y = s.x(ext,s.lenx);
% end

%First check if index is integer or not
rm = t - floor(t); %NOTE don't name this 'rem'! That's a func and will slow things!
if rm == 0
   %integer, so valid index
   y = s.x(ext,t);
else
    %We're between indices, so linear interpolation
    y = s.x(ext,t-rm) * (1-rm) + s.x(ext,t+1-rm) * rm;
    %(t-rm) & (t+1-rm) are much faster than floor(t) & ceil(t)
end
