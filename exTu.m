% This is the example-Tutorial m-file

import gpao.*

maxIter = 50; % number of active-learning iterations

%% Check if we have access to GPML
if ~exist('gp', 'file')
    error('Add GPML in the path please!');
else
    fprintf('Using GPML installed in [%s]\n', which('gp'));
end

%% this is the function we want find the minimum of
% domain represents the domain of the function
[f, domain, trueMinLoc] = testFunctionFactory('f11');
d = size(domain, 1); % dimension of the space

%% initialize the prior
gps = covarianceKernelFactory(1, d);

%% sample a few samples from the Latin Hypercube design
nInit = 7 * d;
obsX = lhsdesign(d, nInit)';
obsY = zeros(size(obsX, 1), 1);
for k = 1:size(obsX, 1)
    obsY(k) = f(obsX(k, :));
end

%% do a litle active learning dance
for k = 1:maxIter
    % ask where to sample next (choose your favorite algorithm)
    %nextX = aoMockus(domain, obsX, obsY, gps);
    nextX = aoKushner(domain, obsX, obsY, gps);

    % evaluate at the suggested point
    nextY = f(nextX);

    % save the measurement pair
    obsX = [obsX; nextX];
    obsY = [obsY; nextY];
end

%% report what has been found
[mv, mloc] = min(obsY);
fprintf('Minimum value: %f found at:\n', mv);
disp(obsX(mloc, :));
fprintf('True minimum value: %f at:\n', f(trueMinLoc));
disp(trueMinLoc)
