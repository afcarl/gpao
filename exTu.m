% This is the example-Tutorial m-file

import gpao.*

maxIter = 50; % number of active-learning iterations
isPlot = true;

%% Check if we have access to GPML
if ~exist('gp', 'file')
    error('Add GPML in the path please!');
else
    fprintf('Using GPML installed in [%s]\n', which('gp'));
end

%% this is the function we want find the minimum of
% domain represents the domain of the function
[f, domain, trueMinLoc] = testFunctionFactory('e1');
d = size(domain, 1); % dimension of the space

%%
if isPlot
    nSeg = 200;
    switch d
        case 1
            x_grid = linspace(domain(1,1), domain(1,2), nSeg)'; % input range
            f_true = f(x_grid);
            fig = figure(5481); clf; subplot(1,2,1); hold all;
            plot(x_grid, f_true);
            plot(trueMinLoc, f(trueMinLoc), 'ro')
            xlabel('domain');
            ylabel('function value');
            subplot(1,2,2);
            xlabel('sample #');
            ylabel('evaluated response');
            line([1, maxIter], f(trueMinLoc) * [1,1], 'LineStyle', ':', 'Color', 'r');
        case 2
            x_grid = linspace(domain(1,1), domain(1,2), nSeg); % input range
            y_grid = linspace(domain(2,1), domain(2,2), nSeg); % input range
            [X, Y] = meshgrid(x_grid, y_grid);
            f_true = f([X(:), Y(:)]);
            f_true = reshape(f_true, nSeg, nSeg);
            nContourLines = 15;
            fig = figure(5481); clf; subplot(1,2,1); hold all
            contour(X, Y, f_true, nContourLines);
        otherwise
            fprintf('Sorry, no support for d > 2\n');
    end
end

%% initialize the prior
gps = covarianceKernelFactory(1, d);

%% sample a few samples from the Latin Hypercube design
nInit = 3 * d;
obsX = lhsdesign(d, nInit)';
o = ones(nInit,1);
obsX = obsX .* (o * (domain(:,2) - domain(:,1))') + o * domain(:,1)';
obsY = zeros(size(obsX, 1), 1);
for k = 1:size(obsX, 1)
    obsY(k) = f(obsX(k, :));
end

if isPlot
    figure(fig); subplot(1,2,1);
    if d == 1
        plot(obsX, obsY, 'bx');
        [ym, ys2, m, s2] = gp(gps.hyp, @infExact, gps.meanfunc, gps.covfunc, ...
            gps.likfunc, obsX, obsY, x_grid);
        
        ph = zeros(3,1);
        ph(1) = plot(x_grid, m, 'r-');
        ph(2) = plot(x_grid, m - 2*s2, 'r:');
        ph(3) = plot(x_grid, m + 2*s2, 'r:');
    else
        plot(obsX(:,1), obsX(:,2), 'bx');
    end
end
               
%% do a litle active learning dance
for k = 1:maxIter
    % ask where to sample next (choose your favorite algorithm)
    %[nextX, gps] = aoMockus(domain, obsX, obsY, gps);
    [nextX, gps] = aoKushner(domain, obsX, obsY, gps);
    
    % evaluate at the suggested point
    nextY = f(nextX) + 0.01 * randn;
    
    % save the measurement pair
    obsX = [obsX; nextX];
    obsY = [obsY; nextY];
    
    %% Plot progress
    if isPlot
        switch d
            case 1
                [ym, ys2, m, s2] = gp(gps.hyp, @infExact, gps.meanfunc, gps.covfunc, ...
                    gps.likfunc, obsX, obsY, x_grid);
                
                figure(fig); subplot(1,2,1);
                plot(nextX, nextY, 'bx');
                set(ph, 'Color', 0.7 * [1,1,1])
                ph(1) = plot(x_grid, m, 'r-');
                ph(2) = plot(x_grid, m - 2*s2, 'r:');
                ph(3) = plot(x_grid, m + 2*s2, 'r:');
            case 2
                subplot(1,2,1);
                plot(nextX(1), nextX(2), 'bx');
        end
    end
    
    if isPlot
        subplot(1,2,2); hold all;
        plot(obsY, 'ko-');
        pause
    end
end

%% report what has been found
[mv, mloc] = min(obsY);
fprintf('Minimum value: %f found at:\n', mv);
disp(obsX(mloc, :));
fprintf('True minimum value: %f at:\n', f(trueMinLoc));
disp(trueMinLoc)
