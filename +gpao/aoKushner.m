function [nextX, gps, xTest, m, s, z, ef] = aoKushner(xRange, observedX, observedY, gps, e)
% [nextX, gps, xTest, m, s, z, ef] = aoKushner(xRange, observedX, observedY, gps)
% Active Optimization Algorithm based on Kushner's method
% Uses GPML library for GP (Gaussian process).
%
% Input:
%   xRange: (d x 2), rectangular Euclidean region
%   observedX: (d x N) N observed x
%   observedY: (N x 1) corresponding y
%   gps: GP structure from previous step (used as initial state for the update)
%	gps.hyp: hyper parameter structure
%	gps.meanfunc: mean function
%	gps.covfunc: covariance function
%	gps.likfunc: likelihood function
%
% Output:
%   nextX: (d x 1) next X to sample from
%   gps: GP structure
%   xTest: set of points where predictive distributions are computed
%   m: mean of predictive distribution
%   s: std of predictive distribution
%   z: current estimate of minimum (usually underestimates)
%   ef: probability density of estimated posterior of Pmin
%   h: entropy of last Pmin
%
% Example parameters: (see GPML documentation by Rasmussen et al.)
% gps.meanfunc = {@meanSum, {@meanLinear, @meanConst}}; gps.hyp.mean = [0; 0];
% gps.covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1;
% gps.hyp.cov = log([ell; sf]);
% gps.likfunc = @likGauss; sn = 0.1; gps.hyp.lik = log(sn);
%
% Algorithm: < Selection of x_{n+1} >
% 1. From GP, estimate z = min(f) via minimizer of mean posterior
%
% Requires: GPML
% See also: runAOsample1.m
%
% Copyright Memming 2011. All rights reserved.
% $Id: aoKushner.m 939 2011-09-22 23:55:08Z memming $

import gpao.*

%% process input
d = size(xRange, 1); % input space dimension
N = size(observedX, 1);
assert(d == size(observedX, 2));
assert(N == numel(observedY)); observedY = observedY(:);

%% Make the grid
nTestPerDim = max(ceil(200^(1/d)), 10); % Make the grid for Pmin sampling
nTest = nTestPerDim^d;

if d == 1 % Grid to evaluate things (1D)
    xTest = linspace(xRange(:,1), xRange(:,2), nTest)'; 
else
    xtr = cell(d, 1);
    for kD = 1:d
	xtr{kD} = linspace(xRange(kD, 1), xRange(kD, 2), nTestPerDim)';
    end
    xTest = zeros(nTest, d);
    for kD = 1:d
	xr = xtr{kD};
	xr = xr(:, ones(1, nTest/nTestPerDim));
	xr = permute(xr, [2:kD 1 kD+1:d]);
	xTest(:, kD) = xr(:);
    end
end

x = observedX; y = observedY;

% Bayesian model selection to find hyperparameter 
gps = evidenceOptHyp(gps, x, y);

%% Perform GP on the test grid
[~, ~, m, s2] = gp(gps.hyp, @infExact, gps.meanfunc, gps.covfunc, ...
		   gps.likfunc, x, y, xTest);
s = sqrt(s2);

%% Estimate the current min f
[z, zi] = min(m);

%% Kushner criterion
% x(next) = arg max_{x} Pr[ f(x) < z - e ]

% z <- z - e;
z = z - s(zi)/10;
% tail probability
ef = erf((z-m)/sqrt(2)./s)/2 + .5; % se = sum(ef); ef = ef / se;

% argmax
[~, idxNext] = max(ef);

nextX = xTest(idxNext, :);

end
