function [rstar, psistar, stability, stabtype] = ...
  rstarfull(a, b1, b2, e, F, W, All)
% [rstar, psistar, stability, stabtype] = 
%     rstarfull(alpha, beta1, beta2, epsilon, force, Omega, All)
%
% Finds r*, psi*, stability (1 or 0) and stability type numerically
% for the full canonical model. Set the optional input argument 'All'
% to 1 (or any nonzero value) to get both stable and unstable fixed points.
% (Default for All is 0, that is, rstarfull outputs only stable points.)
% Stabtype: 4 = stable nodes, 3 = stable spirals, 2 = unstable nodes,
% 1 = unstable spirals, 0 = saddle points

if F < 0
  error('F cannot be negative')
elseif F == 0
  rstar = spontAmp(a,b1,b2,e);
  psistar = zeros(size(rstar));
  stability = ones(size(rstar));
  return
end
if nargin < 7
  All = 0; % default: only stable points
end
if b2 == 0 && e ~= 0
  e = 0;
  %disp('rstarfull: epsilon set to 0 since b2 = 0')
end

% Find roots of steady-state equation numerically
r = sqrt(roots([ (b1-b2)^2*e^2,...
                  -2*(b1-b2)*(b1-a*e)*e,...
                  b1^2-4*a*b1*e+(e*a^2+2*a*b2+e*W^2)*e,...
                  -2*a^2*e+2*a*b1-(e*F^2+2*W^2)*e,...
                  a^2+2*e*F^2+W^2,...
                  -F^2]));
r = r(find(abs(imag(r)) < eps('single'))); % take only real roots
r = real(r);
r = sort(unique(r),'descend'); % remove multiple roots
if b2
  r = r(find(r < 1/sqrt(e))); % take r's below the asymptote
end

% Find corresponding psi's
signPsi = (W >= 0)*2 - 1;
psi = signPsi*acos(-(a*r + b1*r.^3 + b2*e*r.^5./(1-e*r.^2))/F);
maximag = max(abs(imag(psi)));
if  maximag > eps('single')*100
  disp(['Warning (rstarfull): significant nonzero imaginary part in psi ('...
    num2str(maximag) ') for W = ' num2str(W)])
end     
psi = real(psi);

% Jacobian Matrix
J11 = a + 3*b1*r.^2 + e*b2*r.^4.*(5-3*e*r.^2)./(1-e*r.^2).^2;
J12 = -F*sin(psi);
J21 = F*sin(psi)./r.^2;
J22 = -F*cos(psi)./r;
delta = J11.*J22 - J12.*J21; % determinant of Jacobian
tau = J11 + J22; % trace of Jacobian
chdet = tau.^2 - 4*delta; % determinant of characteristic eq

% Stability type
% (4 = stable nodes; 3 = stable spirals; 2 = unstable nodes;
%  1 = unstable spirals; 0 = saddle points)
stabtype = 2*sign(delta) -1*sign(tau) + .5*sign(chdet) + .5;
stabtype(find(stabtype < 0)) = 0;
stability = (stabtype >= 3); % 1 = stable, 0 = unstable

% Output (either stable fixed pts only or both stable and unstable)
if All
  rstar = r;
  psistar = psi;
else
  rstar = r(find(stability));
  psistar = psi(find(stability));
  stability = stability(find(stability));
  stabtype = stabtype(find(stability));
end
