function [nextX, gps, xTest, m, s, z, ef] = aoMockus(xRange, observedX, observedY, gps)
% [nextX, gps, xTest, m, s, z, ef] = aoMockus(xRange, observedX, observedY, gps)
% Active Optimization Algorithm based on Mockus's method
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
%
% Copyright Memming 2011. All rights reserved.
% $Id: aoMockus.m 939 2011-09-22 23:55:08Z memming $

import gpao.*

%% process input
d = size(xRange, 1); % input space dimension
N = size(observedX, 1);
assert(d == size(observedX, 2));
assert(N == numel(observedY)); observedY = observedY(:);

%% Parameters of the algorithm
nMinGridPoints = 200;
isSmoothGrid = true;

%% Make the grid for Pmin sampling
[xTest, xTestDiff, nTest, nTestPerDim] = makeGrid(xRange, nMinGridPoints);

nCandidateSample = nTest; % candidate x_{n+1} to be drawn at every step

x = observedX; y = observedY;

% Bayesian model selection to find hyperparameter 
gps = evidenceOptHyp(gps, x, y);

%% Perform GP on the test grid
[ym ys2 m s2] = gp(gps.hyp, @infExact, gps.meanfunc, gps.covfunc, ...
		   gps.likfunc, x, y, xTest);
s = sqrt(s2);

%% Estimate the current min f
[z, zi] = min(m);

z = z - s(zi)/10;
ef = normcdf2(z, m, s);
%ef = (z-m) .* (erf((z-m)/sqrt(2)./s)/2 + .5) + s .* exp(-(z-m).^2./s2); % se = sum(ef); ef = ef / se; % <== numerically unstable
ef = ef + s .* exp(-(z-m).^2./s2); % se = sum(ef); ef = ef / se;
[dummy, mefi] = max(ef);
idxNext = mefi;

nextX = xTest(idxNext, :);
