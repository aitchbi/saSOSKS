function [cOrds,opts,kType] = hb_get_precomp_cheby_ords(N,kType)
if ~exist('kType','var') || isempty(kType)
    kType = 'MeyerLike'; % type of uniform kernels for estimating energy
end
switch kType
    case 'MeyerLike'
        t = 'meyer';
    case 'SplineType' % not supported yet
        t = 'spline';
end
d = sprintf('cOrds_%s_%d.mat',t,N);
d = which(d);
if ~isempty(d)
    d = load(d);
    cOrds = d.cOrds;
    opts = struct;
    opts.justCheckKernels = false;
else
    cOrds = [];
    opts = struct;
    opts.tol.kernel     = 1e-3;
    opts.tol.tightframe = 1e-3;
    opts.justCheckKernels = false; % Note
    opts.minOrds = 50;
    opts.maxOrds = 800;
    opts.parallelize = true;
end
end

% Note: unif SOSKS should---to a reasonable degree as determined by allowed
% error tolerance specified here, opts_cOrds.tol---form a tight (&
% Parseval) frame for correct estimation of energy spectral density