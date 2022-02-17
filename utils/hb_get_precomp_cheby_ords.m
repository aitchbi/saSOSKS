function [cOrds,opts_cOrds] = hb_get_precomp_cheby_ords(type,N)
switch type
    case 'MeyerLike'
        t = 'meyer';
    case 'SplineType'
        t = 'spline';
end
d = sprintf('cOrds_%s_%d.mat',t,N);
d = which(d);
if ~isempty(d)
    d = load(d);
    cOrds = d.cOrds;
    opts_cOrds = struct;
    opts_cOrds.justCheckKernels = false;
else
    cOrds = [];
    opts_cOrds = struct;
    opts_cOrds.tol.kernel     = 1e-3;
    opts_cOrds.tol.tightframe = 1e-3;
    opts_cOrds.justCheckKernels = false; % Note
    opts_cOrds.minOrds = 50;
    opts_cOrds.maxOrds = 600;
    opts_cOrds.parallelize = true;
end
end

% Note: unif SOSKS should---to a reasonable degree as determined by allowed
% error tolerance specified here---form a tight (& Parseval) frame for
% correct estimation of energy spectral density