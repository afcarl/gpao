function [gps fv1] = evidenceOptHyp(gps, x, y)
% Evidence optimization routine for hyperparamter selection.
% Try previous hyper-parameter as well as a grid of hyper-parameters
% Input
%   gps: GP structure from previous step (used as initial state for the update)
%	gps.hyp: hyper parameter structure
%	gps.meanfunc: mean function
%	gps.covfunc: covariance function
%	gps.likfunc: likelihood function
%   x: (N x d) N observed x
%   y: (N x 1) corresponding y to x
%
% Output
%   gps: optimized parameter set
%   fv1: negative log evidence over iterations of conjugate gradient
%
% Later we may want multiple covariances (TODO)
%
% See also: aoSample2.m
%
% $Id: evidenceOptHyp.m 941 2011-09-23 05:28:48Z memming $
% Copyright Memming 2011. All rights reserved.

mnlh = nan(length(gps.hypgrid), 1); % save marginal negative log likelihood
for k = 1:length(gps.hypgrid)
    mnlh(k) = gp(gps.hypgrid(k), @infExact, gps.meanfunc, gps.covfunc, ...
		gps.likfunc, x, y);
end
[minMlh, idx] = min(mnlh);

% Test previous hyperparameter
fv1 = gp(gps.hyp, @infExact, gps.meanfunc, gps.covfunc, gps.likfunc, x, y);
if fv1 >= minMlh
    fprintf('Grid point is better than previous hyperparameter\n');
    hyp = gps.hypgrid(idx);
else
    hyp = gps.hyp;
end

% hyp = gps.hypgrid(idx);
[hyp1 fv1] = minimize(hyp, @gp, -100, @infExact, gps.meanfunc, ...
		    gps.covfunc, gps.likfunc, x, y);
gps.hyp = hyp1;
