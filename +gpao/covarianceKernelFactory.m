function [gps, desc] = covarianceKernelFactory(setting, d)
% Set of commonly used covariance kernels and their hyperparameter grids
%
% Input
%   setting: numerical selection of options
%   d: input dimension

import gpao.*

maxOption = 7;
if nargin < 1
    gps = maxOption;
    return
end

gps.meanfunc = @meanConst;
gps.hyp.mean = [0];
gps.hypgridRange.mean = -1:1:1;

switch setting
case 1
    innerCovFunc = {@covSEiso};
    desc = 'Gaussian (iso)';
case 2
    innerCovFunc = {@covMaterniso, 1};
    desc = 'Matern (iso) 1';
case 3
    innerCovFunc = {@covMaterniso, 3};
    desc = 'Matern (iso) 3';
case 4
    innerCovFunc = {@covMaterniso, 5};
    desc = 'Matern (iso) 5';
case 5
    innerCovFunc = {@covPPiso, 0};
    desc = 'Piecewise polynomial compact 0 diff';
case 6
    innerCovFunc = {@covPPiso, 1};
    desc = 'Piecewise polynomial compact 2 diff';
case 7
    innerCovFunc = {@covPPiso, 2};
    desc = 'Piecewise polynomial compact 4 diff';
end
% innerCovFunc = {@covRQiso}; % 3 params
% innerCovFunc = {@covSEard}; % requires d+1 hyperparams

covfunc = {@covCaching, innerCovFunc};
covCaching('clear');
covCaching('noverbose');

gps.covfunc = covfunc;
ell = 1; sf = 1;
gps.hyp.cov = log([ell; sf]);
gps.hypgridRange.cov = [(-2:1:0); (-1:1:1)];

gps.likfunc = @likGauss;
sn = 0.1;
gps.hyp.lik = log(sn);
gps.hypgridRange.lik = -3:0;

%% Make the hyper-parameter grid (fixed case for 1 mean, 1 lik, 2 cov)
k = 1;
hyp = gps.hyp;
for k1 = 1:size(gps.hypgridRange.mean, 2)
    hyp.mean = gps.hypgridRange.mean(k1);
    for k4 = 1:size(gps.hypgridRange.lik, 2)
        hyp.lik = gps.hypgridRange.lik(k4);
        p = perms(1:size(gps.hypgridRange.cov, 2));
        for k2 = 1:size(p,1)
            for k3 = 1:size(gps.hypgridRange.cov, 1)
                hyp.cov(k3) = gps.hypgridRange.cov(k3, p(k2, 1));
            end
            gps.hypgrid(k) = hyp;
            k = k + 1;
        end
    end
end
