function K = covCaching(innerCov, hyp, x, z, i)
% Caching interface for GPML
% This covariance wrapper function has internal state for caching kernels and
% it derivatives, and can update kernel matrics by a single element. This is
% useful for online learning.
%
% Input
%   innerCov: function handle or cell type covaraince function
%           (e.g. {@covMaterniso, 3} or @covSEiso)
%           it could also be a string to control the cache state and verbosity
%           'verbose', 'noverbose': debugging messages
%           'reset': clear cache and restart
%   hyp, x, z, i: same as the argument for innerCov
%	    (x: (n x d))
%
% For efficient use of kernel matrices, in the following scenarios:
% S1. x, z both not changed
% S2. z completely changed, but x is the same
% S3. x increased by 1 element
% S4. last element of x is replaced
% while hyp does not change
% If hyp or entire x changes, we recompute.
%
% Usage:
% covCaching maintains internal cached state.
%
% CAUTION:
% It uses 'persistent' keyword, which means you cannot use more than one
% simultaneously.
%
% WARNING:
% z must be fixed for the current implementation.
%
% See also: testCovCaching (unit test)
%
% Copyright 2011 Memming. Same license as GPML.
% $Id: covCaching.m 939 2011-09-22 23:55:08Z memming $

persistent v cached

if ischar(innerCov) % Internal state changing command
    switch lower(innerCov)
        case {'reset', 'clear'}
            clear v cached
        case 'verbose'
            v = true;
        case 'noverbose'
            v = false;
        otherwise
            disp('Unknown command');
    end
    return;
end
if isa(innerCov, 'function_handle'); innerCov = {innerCov}; end
if nargin < 2; K = feval(innerCov{:}); return; end

if isempty(v); v = false; end
if isempty(cached)
    defaultHyp = nan(eval(feval(innerCov{:})), 1); % assume row vector
    cached.Kxx.hyp = defaultHyp; cached.Ks.hyp = defaultHyp;
    cached.kss.hyp = defaultHyp; cached.dKsi.hyp = defaultHyp;
    cached.dKi.hyp = defaultHyp; cached.dkssi.hyp = defaultHyp;
end

assert(size(hyp, 2) == 1, 'Hyper-parameter is assumed to be column vector');
nEffectiveArguments = nargin - 1;

if nEffectiveArguments == 2 || (nEffectiveArguments == 3 && isempty(z))
    if any(cached.Kxx.hyp ~= hyp)
        if v; fprintf('Kxx: recomputing\n'); end
        K = feval(innerCov{:}, hyp, x);
        cached.Kxx.hyp = hyp; cached.Kxx.K = K; cached.Kxx.x = x;
        return;
    end
    
    if size(cached.Kxx.x,1) == size(x,1)-1 ... % adding 1 x
            && all(all(x(1:end-1,:) == cached.Kxx.x))
        if v; fprintf('Kxx: one additional x\n'); end
        Kaa = feval(innerCov{:}, hyp, x(end,:));
        Kxa = feval(innerCov{:}, hyp, x(end,:), cached.Kxx.x);
        K = [cached.Kxx.K Kxa'; Kxa Kaa];
        cached.Kxx.K = K; cached.Kxx.x = x;
        return;
    elseif size(cached.Kxx.x,1) == size(x,1)
        if all(all(x == cached.Kxx.x)) % identical x
            if v; fprintf('Kxx: using cached\n'); end
            K = cached.Kxx.K;
            return;
        else
            if size(x,1) > 1 && ...
                    all(all(x(1:end-1,:) == cached.Kxx.x(1:end-1,:))) 
		% last one replaced
                if v; fprintf('Kxx: last x replaced\n'); end
                Kaa = feval(innerCov{:}, hyp, x(end,:));
                Kxa = feval(innerCov{:}, hyp, x(end,:), cached.Kxx.x(1:end-1,:));
                K = [cached.Kxx.K(1:end-1,1:end-1) Kxa'; Kxa Kaa];
                cached.Kxx.K = K; cached.Kxx.x = x;
                return;
            end
        end
    end
    
    % if you reached here, we recompute
    if v; fprintf('Kxx: recomputing\n'); end
    K = feval(innerCov{:}, hyp, x);
    cached.Kxx.K = K; cached.Kxx.x = x;
    return;
end

%% kss = covNAME(hyp, xs, 'diag')
% For constant diagonal (normalized) covariance kernels, doing the caching
% is slower because of the overhead.
if nEffectiveArguments == 3 && ischar(z) && strcmp(z, 'diag')
    if any(cached.kss.hyp ~= hyp)
        K = feval(innerCov{:}, hyp, x, 'diag');
        cached.kss.hyp = hyp; cached.kss.K = K; cached.kss.x = x;
        return;
    end
    
    if size(cached.kss.x,1) == size(x,1)-1 ... % adding 1 x
            && all(all(x(1:end-1,:) == cached.kss.x))
        Kaa = feval(innerCov{:}, hyp, x(end,:), 'diag');
        K = [cached.kss.K Kaa];
        cached.kss.hyp = hyp; cached.kss.K = K; cached.kss.x = x;
        return;
    elseif size(cached.kss.x,1) == size(x,1)
        if all(all(x == cached.kss.x)) % identical x
            K = cached.kss.K;
            return;
        else
            if size(x,1) > 1 && ...
                    all(all(x(1:end-1,:) == cached.kss.x(1:end-1,:))) % last one replaced
                Kaa = feval(innerCov{:}, hyp, x(end,:), 'diag');
                K = [cached.kss.K(1:end-1) Kaa];
                cached.kss.hyp = hyp; cached.kss.K = K; cached.kss.x = x;
                return;
            end
        end
    end
    % if you reached here, we recompute
    K = feval(innerCov{:}, hyp, x, 'diag');
    cached.kss.hyp = hyp; cached.kss.K = K; cached.kss.x = x;
    return;
end

%% Ks  = covNAME(hyp, x, xs)
if nEffectiveArguments == 3
    if any(cached.Ks.hyp ~= hyp)
        if v; fprintf('Kxz: recomputing\n'); end
        K = feval(innerCov{:}, hyp, x, z);
        cached.Ks.hyp = hyp; cached.Ks.K = K; cached.Ks.x = x; cached.Ks.z = z;
        return;
    end
    
    % TODO: here we assume z is fixed, only x chages
    % assert(all(all(z == cached.Ks.z)), 'z must be the same');
    if any(size(z) ~= size(cached.Ks.z)) || ~all(all(z == cached.Ks.z));
        if v; fprintf('Kxz: recomputing (z change)\n'); end
        K = feval(innerCov{:}, hyp, x, z);
        cached.Ks.hyp = hyp; cached.Ks.K = K; cached.Ks.x = x; cached.Ks.z = z;
        return;
    end
    
    if size(cached.Ks.x,1) == size(x,1)-1 ... % adding 1 x
            && all(all(x(1:end-1,:) == cached.Ks.x))
        if v; fprintf('Kxz: one additional x\n'); end
        Kxz = feval(innerCov{:}, hyp, x(end,:), z);
        K = [cached.Ks.K; Kxz];
        cached.Ks.hyp = hyp; cached.Ks.K = K; cached.Ks.x = x; cached.Ks.z = z;
        return;
    elseif size(cached.Ks.x,1) == size(x,1)
        if all(all(x == cached.Ks.x)) % identical x
            if v; fprintf('Kxz: using cached\n'); end
            K = cached.Ks.K;
            return;
        else
            if size(x,1) > 1 && ...
                    all(all(x(1:end-1,:) == cached.Ks.x(1:end-1,:))) % last one replaced
                if v; fprintf('Kxz: last x replaced\n'); end
                Kxz = feval(innerCov{:}, hyp, x(end,:), z);
                K = [cached.Ks.K(1:end-1,:); Kxz];
                cached.Ks.hyp = hyp; cached.Ks.K = K;
                cached.Ks.x = x; cached.Ks.z = z;
                return;
            end
        end
    end
    
    % if you reached here, we recompute
    if v; fprintf('Kxz: recomputing\n'); end
    K = feval(innerCov{:}, hyp, x, z);
    cached.Ks.hyp = hyp; cached.Ks.K = K; cached.Ks.x = x; cached.Ks.z = z;
    return;
end

%% dKi, dKsi, dkssi
assert(nEffectiveArguments == 4)
if isempty(z) % dKi
    if size(cached.dKi,1) < i || any(cached.dKi(i).hyp ~= hyp)
        if v; fprintf('dKi: recomputing\n'); end
        K = feval(innerCov{:}, hyp, x, [], i);
        cached.dKi(i).hyp = hyp; cached.dKi(i).K = K; cached.dKi(i).x = x;
        return;
    end
    
    if size(cached.dKi(i).x,1) == size(x,1)-1 ... % adding 1 x
            && all(all(x(1:end-1,:) == cached.dKi(i).x))
        if v; fprintf('dKi: one additional x\n'); end
        Kaa = feval(innerCov{:}, hyp, x(end,:), [], i);
        Kxa = feval(innerCov{:}, hyp, x(end,:), cached.dKi(i).x, i);
        K = [cached.dKi(i).K Kxa'; Kxa Kaa];
        cached.dKi(i).K = K; cached.dKi(i).x = x;
        return;
    elseif size(cached.dKi(i).x,1) == size(x,1)
        if all(all(x == cached.dKi(i).x)) % identical x
            if v; fprintf('dKi: using cached\n'); end
            K = cached.dKi(i).K;
            return;
        else
            if size(x,1) > 1 && ...
                    all(all(x(1:end-1,:) == cached.dKi(i).x(1:end-1,:))) % last one replaced
                if v; fprintf('dKi: last x replaced\n'); end
                Kaa = feval(innerCov{:}, hyp, x(end,:), [], i);
                Kxa = feval(innerCov{:}, hyp, x(end,:), ...
                    cached.dKi(i).x(1:end-1,:), i);
                K = [cached.dKi(i).K(1:end-1,1:end-1) Kxa'; Kxa Kaa];
                cached.dKi(i).K = K; cached.dKi(i).x = x;
                return;
            end
        end
    end
    
    % if you reached here, we recompute
    if v; fprintf('dKi: recomputing\n'); end
    K = feval(innerCov{:}, hyp, x, [], i);
    cached.dKi(i).K = K; cached.dKi(i).x = x;
    return;
elseif ischar(z) && strcmp(z, 'diag') % dkssi
    if size(cached.dkssi,1) < i || any(cached.dkssi(i).hyp ~= hyp)
        if v; fprintf('dkssi: recomputing (new/hyp)\n'); end
        K = feval(innerCov{:}, hyp, x, 'diag', i);
        cached.dkssi(i).hyp = hyp; cached.dkssi(i).K = K; cached.dkssi(i).x = x;
        return;
    end
    
    % TODO: here we assume z is fixed, only x chages
    % note that x is z in this case for GP
    assert(all(all(x == cached.dkssi(i).x)), 'z change not implemented');
    if v; fprintf('dkssi: recomputing\n'); end
    K = cached.dkssi(i).K;
    return;
else % dKsi
    if size(cached.dKsi,1) < i || any(cached.dKsi(i).hyp ~= hyp)
        if v; fprintf('dKsi: recomputing (new/hyp)\n'); end
        K = feval(innerCov{:}, hyp, x, z, i);
        cached.dKsi(i).hyp = hyp; cached.dKsi(i).K = K;
        cached.dKsi(i).x = x; cached.dKsi(i).z = z;
        return;
    end
    
    % TODO: here we assume z is fixed, only x chages
    % assert(all(all(z == cached.dKsi(i).z)), 'z must be the same');
    if any(size(z) ~= size(cached.Ks.z)) || ~all(all(z == cached.dKsi(i).z));
        if v; fprintf('dKsi: recomputing (new/hyp)\n'); end
        K = feval(innerCov{:}, hyp, x, z, i);
        cached.dKsi(i).hyp = hyp; cached.dKsi(i).K = K;
        cached.dKsi(i).x = x; cached.dKsi(i).z = z;
        return;
    end
    
    if size(cached.dKsi(i).x,1) == size(x,1)-1 ... % adding 1 x
            && all(all(x(1:end-1,:) == cached.dKsi(i).x))
        if v; fprintf('dKsi: adding 1 x\n'); end
        Kxz = feval(innerCov{:}, hyp, x(end,:), z, i);
        K = [cached.dKsi(i).K; Kxz];
        cached.dKsi(i).K = K; cached.dKsi(i).x = x; cached.dKsi(i).z = z;
        return;
    elseif size(cached.dKsi(i).x,1) == size(x,1)
        if all(all(x == cached.dKsi(i).x)) % identical x
            if v; fprintf('dKsi: cache exact match\n'); end
            K = cached.dKsi(i).K;
            return;
        else
            if size(x,1) > 1 && ...
                    all(all(x(1:end-1,:) == cached.dKsi(i).x(1:end-1,:))) % last one replaced
                if v; fprintf('dKsi: last 1 x replaced\n'); end
                Kxz = feval(innerCov{:}, hyp, x(end,:), z, i);
                K = [cached.dKsi(i).K(1:end-1,:); Kxz];
                cached.dKsi(i).K = K; cached.dKsi(i).x = x; cached.dKsi(i).z =z;
                return;
            end
        end
    end
    
    % if you reached here, we recompute
    if v; fprintf('dKsi: recomputing\n'); end
    K = feval(innerCov{:}, hyp, x, z, i);
    cached.dKsi(i).hyp = hyp; cached.dKsi(i).K = K;
    cached.dKsi(i).x = x; cached.dKsi(i).z = z;
    return;
end

error('Memming''s wrong');
