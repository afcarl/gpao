import gpao.*

functionNameList = {'e1', 'e2', 'e3', 'e4', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11'};
%functionNameList = {'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11'}; % 2D functions only
%functionNameList = {'e1', 'e2', 'e3', 'e4'}; % 1D functions only

for kFunctionName = 1:length(functionNameList)
    functionName = functionNameList{kFunctionName};
    disp(['>>> ' functionName ' <<<']);
    [f, domain, trueMinLoc] = testFunctionFactory(functionName);
    d = size(domain, 1);

    clf;
    if d == 1
	x = linspace(domain(1,1), domain(1,2))';
	hold on;
	plot(x, f(x), '-');
	for kTrueMin = 1:size(trueMinLoc,2)
	    tx = trueMinLoc(kTrueMin);
	    fv = f(tx);
	    plot(tx, fv, 'r.-');
	end
	if exist('fminbnd', 'file')
	    [xhat, fv] = fminbnd(f, domain(1), domain(2));
	    plot(xhat, fv, 'kd', 'MarkerSize', 10);
	end
    elseif d == 2
	x = linspace(domain(1,1), domain(1,2), 151);
	y = linspace(domain(2,1), domain(2,2), 151);
	[X, Y] = meshgrid(x, y);
	fXY = f([X(:), Y(:)]);
	fXY = reshape(fXY, 151, 151);
	subplot(1,2,2); cla; hold on;
	mesh(X, Y, fXY);
	
	subplot(1,2,1); cla; hold on;
	contour(X, Y, fXY, 25);
	empMin = min(fXY(:));
	for kTrueMin = 1:size(trueMinLoc,2)
	    tx = trueMinLoc(1, kTrueMin); ty = trueMinLoc(2, kTrueMin);
	    fv = f([tx, ty]);
	    fprintf('At (%g, %g), value (%g). Empirical global min [%g]\n', ...
		tx, ty, fv, empMin);
	    if fv > empMin
		fprintf('Something''s wrong\n');
	    end
	    subplot(1,2,1);
	    plot3(tx, ty, fv, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0]);
	    subplot(1,2,2);
	    plot3(tx, ty, fv, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0]);
	end
	axis tight
    else
	error('Higher dimension visualization is not implemented');
    end
    title(functionName);
    disp('Press any key...');
    pause;
end
