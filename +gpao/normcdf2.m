function [p] = normcdf2(x, mu, sigma)
% p = normcdf2(x, mu, sigma)
% CDF of Gaussian distribution, often denoted as $\Phi$.
% normcdf in stat toolbox has more functionality.
%
% Input
%   x: (N x 1) points to evaluate normcdf
%   mu: (N x 1) mean
%   sigma: (N x 1) standard deviation
%
% Output
%   p: (N x 1) probability
%
% $Id: normcdf2.m 938 2011-09-22 20:05:40Z memming $
% Copyright Memming 2011. All rights reserved.

% erfc is better than erf for this purpose
% adding large number to small number is numerically unstable!!!
% p = erf((x-m)/sqrt(2)./s)/2 + .5;
z = (x - mu) ./ sigma;
p = erfc(-z ./ sqrt(2)) / 2;
