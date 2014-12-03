%% spontAmp
%   r = spontAmp(a, b1, b2, e)
%
%  Finds spontaneous amplitude(s) of fully expanded canonical model
%%
function r = spontAmp(a, b1, b2, e)

if b2 == 0 && e ~=0
  e = 0;
  %disp('spontAmp: epsilon set to 0 since b2 = 0')
end

% Find r* numerically
r = roots([e*(b2-b1), 0, b1-e*a, 0, a, 0]);
r = real(unique(r(find(abs(imag(r)) < eps('single')))));
                    % only unique real values
r = r(find(r >=0)); % no negative amplitude
if b2
  r = r(find(r < 1/sqrt(e))); % r* below the asymptote only
end

% Take only stable r*
ind1 = find(slope(r,a,b1,b2,e) < 0);
ind2a = find(slope(r,a,b1,b2,e) == 0);
ind2b = find(slope(r-eps('single'),a,b1,b2,e) < 0);
ind2c = find(slope(r+eps('single'),a,b1,b2,e) < 0);
ind2 = intersect(ind2a,intersect(ind2b,ind2c));
r = r([ind1; ind2]);
r = sort(r,'descend');

% ========================================================
function drdotdr = slope(r, a, b1, b2, e)
drdotdr = a + 3*b1*r.^2 + (5*e*b2*r.^4-3*e^2*b2*r.^6)./((1-e*r.^2).^2);