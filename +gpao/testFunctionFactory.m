function [f, domain, trueMinLoc] = testFunctionFactory(name)
% Returns an oracle (without noise)
% This is not a unit test of another function, but rather a factory for
% "test functions".
%
% Output
%   f: @(x) oracle function handle
%   domain: (d x 2) range of optimization (min and max for each dimension)
%   trueMinLoc: (d x n) n true known global minima 
%		(could be analytical or numerically optained)
%
% Reference: Cox & John
% Reference: http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page364.htm
% $Id: testFunctionFactory.m 939 2011-09-22 23:55:08Z memming $
% Copyright 2011 Memming. All rights reserved.

switch lower(name)
case {'e1', 'sin-sqrt'} % 1 local minima
    f = @(x) sin(x * 2 * pi * 2.5) .* sqrt(x);
    domain = [0 1];
    trueMinLoc = 0.702878196278343;
case {'e2', 'sin'}
    f = @(x) sin(x * 2 * pi * 2.5);
    domain = [0 1];
    trueMinLoc = [0.3, 0.7];
case {'e3'} % 1 local minima
    f = @(x) -(x.^2 + 1) .* cos(3*pi*x);
    domain = [-1 1];
    % wolframalpha query "minimize -(x^2 + 1) * cos(3*pi*x)"
    trueMinLoc = 0.677086 * [-1 1];
case {'e4'} % 4 local minima
    f = @(x) (-exp(-x.^2) + 1) .* cos(3*pi*x);
    domain = [-1.5 1.5];
    % wolframalpha query 
    %	"minimize (-exp(-x^2) + 1) * cos(3*pi*x) from x=-1.5 to 1.5"
    trueMinLoc = 1.01269 * [-1 1];
case {'f1', 'simulated annealing i'}
    % Cox & John f1
    f = @(x,y)(2*x.^2 + 2*y.^2 - 0.3*cos(3*pi*x) - 0.4*cos(4*pi*y) + 0.7);
    domain = [-1 1; -1 1];
    trueMinLoc = [0; 0];
case {'f2', 'simulated annealing ii'}
    % Cox & John f2
    f = @(x,y)(2*x.^2 + 2*y.^2 - 0.3*cos(3*pi*x).*cos(4*pi*y) + 0.3);
    domain = [-1 1; -1 1];
    trueMinLoc = [0; 0];
case {'f3', 'simulated annealing iii'}
    % Cox & John f3
    f = @(x,y)(2*x.^2 + 2*y.^2 - 0.3*cos(3*pi*x + 4*pi*y) + 0.3);
    domain = [-1 1; -1 1];
    trueMinLoc = [0; 0];
case {'f4', 'sine function'}
    % Cox & John f4
    f = @(x,y)(sin(x).^2 + sin(y).^2 - 0.1 * exp(-x.^2 - y.^2));
    domain = [-10 10; -10 10];
    trueMinLoc = [0; 0];
case {'f5', 'tree hump camel-back'}
    f = @(x,y)(2*x.^2 - 1.05 * x.^4 + x.^6 / 6 - x.*y + y.^2);
    domain = [-3 3;-1.5 1.5];
    trueMinLoc = [0; 0];
case {'f6', 'hosaki'}
    fx = @(x)(1 - 8*x + 7*x.^2 - (7/3)*(x.^3) + (x.^4)/4);
    fy = @(y)((y.^2) .* exp(-y));
    f = @(x,y) fx(x) .* fy(y);
    domain = [0 5;0 6];
    trueMinLoc = [4; 2];
    localMinLoc = [1; 2];
case {'f7', 'goldstein and price'}
%    f = @(x,y)((1+(x+y+1).^2 ...
%	    * (19 - 14*x + 3*x.^2 - 14*y + 6*x.*y + 3*y.^2)) ...
%	* (30 + (2*x - 3*y).^2 ...
%	    .* (18 - 32*x + 12*x.^2 + 48*y - 36*x.*y + 27*y.^2)));
    % Lizotte Ph.D. dissertation p. 82
    a = @(x,y) 1+((x+y+1).^2).*(19-14*x+3*x.^2-14*y+6*x.*y+3*y.^2);
    b = @(x,y) 30+((2*x-3*y).^2).*(18-32*x+12*x.^2+48*y-36*x.*y+27*y.^2);
    f = @(x,y) a(x,y).*b(x,y);
    % Schonlau Ph.D. dissertation p. 51
    domain = [-2 2;-2 2];
    trueMinLoc = [0; -1];
case {'f8', 'branin'}
    %f = @(x,y)((y - (5/4)*(pi^2)*x.^2 + (5/pi)*x - 6).^2 + 10 * (1 - 1/(8*pi)) * cos(x) + 10);
    % Ref: Torn and Zilinskas 1989
    f = @(x,y)(y-(5.1/(4*pi^2))*x.^2+5*x/pi-6).^2+10*(1-1/(8*pi))*cos(x)+10;
    domain = [-5 10; 0 15];
    trueMinLoc = [-pi, pi, 9.4248; 12.25, 2.25, 2.475]';
case {'f9', 'dixon-szego'}
    % Ref: http://www2.maths.ox.ac.uk/chebfun/examples/opt/html/DixonSzego.shtml
    % L.C.W. Dixon and G.P. Szegő, The global optimization problem: an introduction, in L.C.W. Dixon and G.P. Szegő (eds.), Towards Global Optimisation 2, North-Holland, Amsterdam 1978, pp. 1-15.
    f = @(x,y) (4-2.1*x.^2+ x.^4/3).*x.^2 + x.*y + 4*(y.^2-1).*y.^2;
    domain = [-2 2; -1.25 1.25];
    trueMinLoc = [0.0898, -0.0898; -0.7127, 0.7127]'; % numerical result
case {'f10', 'rosenbrock'}
    % H. H. Rosenbrock, "An automatic method for finding the greatest or least value of a function", Computer Journal 3 (1960), 175-184.
    f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
    domain = [-1.5 1.5; -1 3]';
    trueMinLoc = [1; 1];
case {'f11', 'peaks'}
    % MATLAB Peaks function
    f = @(x,y) peaks(x,y);
    domain = [-3 3; -3 3];
    trueMinLoc = [0.25; -1.625];
otherwise
    error('Unknown test function name [%s]', name);
end

d = size(domain, 1);
if d == 1
    f = @(x)(f(x(:,1)));
elseif d == 2
    f = @(x)(f(x(:,1), x(:,2)));
else
    error('High dimension?');
end
