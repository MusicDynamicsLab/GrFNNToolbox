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
function [dzdt] = zdotAdaptive(t, zColumn, M)

% numberNetworks = length(M.n);
% zMat = cell(numberNetworks,1);
% dzdtMat = cell(numberNetworks,1);
% ind = 1;
% zColumnLength = 0;
% for i = 1:numberNetworks
%     zMat{i} = zColumn(ind:M.n{i}.N + ind -1);
%     ind = ind + M.n{i}.N;
%     zColumnLength = zColumnLength + M.n{i}.N;
% end

nOsc = length(zColumn)/2;

%% Initialize variables and stimulus

% for nx = 1:numberNetworks
    
    z = zColumn(1:nOsc);
    
    n   = M.n{1};
    % z   = n.z;
    a   = n.a;
    b1  = n.b1;
    b2  = n.b2;
    e   = n.e;
    
    %% Add input
    
    x = 0;
    
    for cx = 1:length(n.con)
        con = n.con{cx};
        y = stimulusRun(M.s{1}.x,t*M.fs+1);
        x = x + con.w.*(con.C*y);
    end
    
    %% The differential equation
    % $\dot{z} = z \left( \alpha + \textrm{i}\omega + (\beta_1 + \textrm{i}\delta_1) |z|^2 + \frac{\epsilon (\beta_2 + \textrm{i}\delta_2) |z|^4}{1-\epsilon |z|^2} \right) + x$
    dzdt1 = z.*(a + b1.*abs(z).^2 + e*b2.*(abs(z).^4)./(1-e*abs(z).^2)) + x;
    
% end


%%




%% Initialize variables and stimulus

% for nx = 1:numberNetworks
    
    z = zColumn(nOsc+1:end);
    
    n   = M.n{2};
    % z   = n.z;
    a   = n.a;
    b1  = n.b1;
    b2  = n.b2;
    e   = n.e;
    
    %% Add input
    
    x = 0;
    
    for cx = 1:length(n.con)
        con = n.con{cx};
       
        y = zColumn(1:nOsc);
        
%                 if con.no11
%                     x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ...
%                         - A(e^2, z*y').*repmat(y.'.*A(e^2, abs(y.').^2), con.targetN, 1) ), 2);
%                 else
% %                     x =  x + con.w .* sum(con.C.*( A(e, z)*P_new(e, y.') ), 2);
                    x =  x + con.w .* (con.C*( sum(A(e, z)*P_new(e, y.'),2) ));
%                 end
%                 x = x + con.w .* (con.C*sum(y+sqrt(e)*y.*conj(y)));
%                 x = x + con.w .* (con.C*sum(y));

             
    end
    
    %% The differential equation
    % $\dot{z} = z \left( \alpha + \textrm{i}\omega + (\beta_1 + \textrm{i}\delta_1) |z|^2 + \frac{\epsilon (\beta_2 + \textrm{i}\delta_2) |z|^4}{1-\epsilon |z|^2} \right) + x$
    dzdt2 = z.*(a + b1.*abs(z).^2 + e*b2.*(abs(z).^4)./(1-e*abs(z).^2)) + x;
    dzdt=[dzdt1;dzdt2];
    
% end







% dzdt = NaN(zColumnLength,1);
% ind = 1;
% for i = 1:numberNetworks
%     dzdt(ind:ind + M.n{i}.N - 1) = dzdtMat{i};
%     ind = ind + M.n{i}.N;
% end

%%
function s = stimulusRun(s, time)

% ratio = (tx) - round(tx);
% if floor(tx) < length(s)
%     disp(floor(tx));disp(tx);disp(ceil(tx))
%     s = s(:,floor(tx) + 1) + ratio * (s(:,ceil(tx) + 1) - s(:,floor(tx) + 1));
% else
%     s = s(:,floor(tx) + 1);
% end


% rm = tx - floor(tx); 
% if rm == 0  % integer, so valid index
%     s = s(:,tx);
% else        % We're between indices, so linear interpolation
%     s = s(:,round(tx-rm)) * (1-rm) + s(:,round(tx+1-rm)) * rm;
% end

leftInd = floor(time);
here    = s(:,leftInd);
% disp(time)
% disp(ceil(time))
there   = s(:,ceil(time));
s = here + ((time - leftInd) * (there-here));


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
