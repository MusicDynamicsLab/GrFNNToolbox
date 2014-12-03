%% relTimeDriven11
%   [relTime, r0, rStar] = relTimeDriven11(alpha, beta1, beta2, epsilon, F)
%
%  Calculates an appoximation of the relaxation time (relTime)
%  for a canonical oscillator driven by a sinusoid at its natural frequency
%  (Omega = 0), starting from its spontaneous amplitude (r0) until it
%  reaches within 1 of the driven steady-state amplitude (rStar). When
%  there are multiple spontaneous amplitudes, vectors are given as output.
%  Only dr/dt equation is used, assuming psi = 0 throughout the trajectory.

function [relTime, r0, rStar] = relTimeDriven11(a, b1, b2, e, F)

r0 = spontAmp(a, b1, b2, e); % spontaneous amplitude
rss = rStarDriven11(a, b1, b2, e, F, 0, 0); % stable steady-state amp
rss = sort(rss); % sort in ascending order
rStar = zeros(size(r0));
relTime = zeros(size(r0));

for n = 1:length(r0)
  rssBig = rss(rss > r0(n)); % driven amplitudes larger than r0
  rStar(n) = rssBig(1); % take the one closest to r0
  R = linspace(r0(n),rStar(n),101); % 100 amplitude steps between r0
                                    % and rStar
  dr = R(2)-R(1); % step size
  relTime(n) = dr*sum(dtdr(a, b1, b2, e, F, R(1:100))); % Euler integration
end

function dtdr = dtdr(a, b1, b2, e, F, r)
dtdr = 1./(a*r + b1*r.^3 + e*b2*r.^5./(1-e*r.^2) + F);
