%% rStarDriven11Learn
%   [rStar, AStar, psiStar, stability, stabType] = rStarDriven11Learn(alpha, beta1, lambda, mu1, kappa, forcing amplitude, Omega, All)
%
%  Finds r*, A*, psi*, stability (1 or 0), and stability type (0-4)
%  numerically for an oscillator (truncated after the cubic term) driven by
%  a sinusoidal input via a plastic 1:1 coupling. Omega is the difference
%  between the oscillator's natural frequency and the input frequency
%  in radian.
%  Set the optional input argument 'All' to 1 (or any nonzero value)
%  to get both stable and unstable fixed points. (Default for All is 0,
%  that is, rStarDriven11Learn outputs only stable fixed points.)
%
%  stability: 1 = stable, 0 = unstable
%
%  stabType: 4 = stable node, 3 = stable spiral, 2 = unstable node,
%  1 = unstable spiral, 0 = saddle point

%% Equation
% $$\frac{dz}{dt} = z\left(\alpha + \textrm{i}\omega + \beta_1 |z|^2\right)
% + cx$$
% 
% $$\frac{dc}{dt} = c\left(\lambda + \mu_1 |c|^2\right)
% + \kappa z \bar{x}$$
%
% where $z = re^{\textrm{i}\phi}, c = Ae^{\textrm{i}\theta},
% x = Fe^{\textrm{i}\vartheta}, \psi = \theta - \phi + \vartheta,
% \vartheta = \omega_0t + \theta_0,$ and $\Omega = \omega - \omega_0$

function [rStar, AStar, psiStar, stability, stabType] = ...
  rStarDriven11Learn(a, b1, l, m1, ka, F, W, All)

warning = 0;

%% Check input arguments
if F <= 0
  error('F must be positive')
end
if nargin < 8
  All = 0; % default: only stable points
end

%% Get steady-state connection amplitudes numerically
A = sqrt(roots([-(power(b1,5)*power(ka,3)*power(m1,3)*...
  power(b1 - ka*m1,2)),...
  -(power(b1,5)*power(ka,3)*power(m1,2)*(b1 - ka*m1)*...
  (3*b1*l - ka*(2*a + 5*l)*m1)),...
  -(power(b1,5)*power(ka,3)*m1*(3*power(b1,2)*power(l,2) - ...
  2*b1*ka*m1*(3*a*l + 6*power(l,2) - power(W,2)) + ...
  power(ka,2)*power(m1,2)*(power(a,2) + ...
  8*a*l + 2*(5*power(l,2) + power(W,2))))),...
  -(power(b1,4)*power(ka,3)*(power(b1,3)*power(l,3) + ...
  a*power(F,2)*power(ka,4)*power(m1,3) + ...
  power(b1,2)*ka*m1*(a*power(F,2)*ka - 6*a*power(l,2) - ...
  8*power(l,3) + 4*l*power(W,2)) + b1*power(ka,2)*power(m1,2)*...
  (3*power(a,2)*l + 10*power(l,3) + 6*l*power(W,2) + ...
  2*a*(-(power(F,2)*ka) + 6*power(l,2) + power(W,2))))),...
  power(b1,4)*power(ka,4)*(power(F,2)*power(ka,3)*...
  (-2*power(a,2) + power(F,2)*ka - 3*a*l)*power(m1,2) + ...
  power(b1,2)*(power(F,4)*power(ka,2) - a*power(F,2)*ka*l + ...
  2*power(l,2)*(a*l + power(l,2) - power(W,2))) - ...
  b1*ka*m1*(2*power(F,4)*power(ka,2) + 5*power(l,4) - ...
  4*power(F,2)*ka*power(W,2) + 6*power(l,2)*power(W,2) + ...
  power(W,4) + 4*a*l*(-(power(F,2)*ka) + 2*power(l,2) + ...
  power(W,2)) + power(a,2)*...
  (-2*power(F,2)*ka + 3*power(l,2) + power(W,2)))),...
  -(power(b1,4)*power(ka,5)*(power(a,3)*power(F,2)*power(ka,2)*m1 + ...
  power(a,2)*l*(4*power(F,2)*power(ka,2)*m1 + ...
  b1*(-2*power(F,2)*ka + power(l,2) + power(W,2))) + ...
  l*(-2*power(F,4)*power(ka,3)*m1 + b1*(2*power(F,4)*power(ka,2) - ...
  4*power(F,2)*ka*power(W,2) + power(power(l,2) + power(W,2),2))) + ...
  a*(power(F,2)*power(ka,2)*m1*(-2*power(F,2)*ka + 3*power(l,2) + ...
  power(W,2)) + b1*(2*power(F,4)*power(ka,2) + 2*power(l,2)*(power(l,2)...
  + power(W,2)) - power(F,2)*ka*(2*power(l,2) + power(W,2)))))),...
  power(b1,4)*power(F,2)*power(ka,7)*(-(power(a,3)*l) + ...
  power(F,2)*ka*power(l,2) + power(a,2)*(power(F,2)*ka - 2*power(l,2))...
  -a*l*(-2*power(F,2)*ka + power(l,2) + power(W,2)))]));
A = A(imag(A)==0);

%% Get oscillator amplitudes
r = [sqrt((-ka*a+sqrt(ka^2*a^2+4*ka*b1*A.^2.*(l+m1*A.^2)))/(2*ka*b1));...
  sqrt((-ka*a-sqrt(ka^2*a^2+4*ka*b1*A.^2.*(l+m1*A.^2)))/(2*ka*b1))];
A = [A;A];
ind = find(imag(r)==0);
r = r(ind); A = A(ind);

%% Get relative phases
signPsi = -((W >= 0)*2 - 1);
psi = signPsi*acos(-(a*r+b1*r.^3)./(F*A));

%% Eliminate extraneous solutions
ind = find(abs(a*r+b1*r.^3+F*A.*cos(psi)) < eps('single'));
r = r(ind); A = A(ind); psi = psi(ind);
ind = find(abs(-W-F*(ka*r.^2+A.^2)./(r.*A).*sin(psi)) < eps('single'));
r = r(ind); A = A(ind); psi = psi(ind);

if warning
  maxImag = max(abs(imag(psi)));
  if  maxImag > eps('single')*100
    disp(['Warning (rStarDriven11Learn): significant nonzero'...
      ' imaginary part in psi (' num2str(maxImag) ') for W = ' num2str(W)])
  end
end
psi = real(psi);

%% Get stability type from Jacobian matrix
stabType = zeros(size(r));
for n = 1:length(r)
  J = zeros(3);
  J(1,1) = a + 3*b1*r(n)^2;
  J(1,2) = F*cos(psi(n));
  J(1,3) = -F*A(n)*sin(psi(n));
  J(2,1) = ka*F*cos(psi(n));
  J(2,2) = l + 3*m1*A(n)^2;
  J(2,3) = -ka*F*r(n)*sin(psi(n));
  J(3,1) = -F*(ka/A(n)-A(n)/r(n)^2)*sin(psi(n));
  J(3,2) = F*(ka*r(n)/A(n)^2 - 1/r(n))*sin(psi(n));
  J(3,3) = -F*(ka*r(n)^2 + A(n)^2)/(r(n)*A(n))*cos(psi(n));
  if ~any(isnan(J(:)))
    ev = eig(J);
    if isreal(ev) && all(ev < 0)
      stabType(n) = 4; % stable node
    elseif all(real(ev) < 0)
      stabType(n) = 3; % stable spiral
    elseif isreal(ev) && all(ev > 0)
      stabType(n) = 2; % unstable node
    elseif all(real(ev) > 0)
      stabType(n) = 1; % unstable spiral
    end % saddle pt otherwise
  end
end
stability = (stabType >= 3); % 1 = stable, 0 = unstable

%% Prepare output
if All % both stable and unstable fixed points
  rStar = r;
  AStar = A;
  psiStar = psi;
else % only stable points
  indStab = find(stability);
  rStar = r(indStab);
  AStar = A(indStab);
  psiStar = psi(indStab);
  stability = stability(indStab);
  stabType = stabType(indStab);
end
