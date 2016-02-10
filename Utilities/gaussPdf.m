%% gaussPdf
%  y = gaussPdf(x,mu,sigma)
%
%  Creates a portion y of a Gaussian probability density function along the
%  x axis specified, with mean mu and standard deviation sigma.

function y = gaussPdf(x,mu,sigma)

y=exp(-.5*((x-mu)./sigma).^2)./(sqrt(2*pi).*sigma);