function [coeff,cOrds,G] = hb_graph_filt(G,s,g,cOrds,opts_cOrds)
% HB_GRAPH_FILT filters a given graph signal using the set of spectral graph
% kernels given in g using Chebyshev polynomial approximation. 
%
% Inputs:
%   G : graph structure, at minimal with field A (adjacency matrix)
%   g : cell array of spectral kernels. 
%   cOrds: structure with fields'kernel' and 'tightframe' each of
%   length(g), specifying the Chebyshev polynomial order to use for
%   appriximating each kernel in g. If [], polynomial orders will be
%   estimated.
%
% Outputs: 
%   coeff: filtered signal.
%   cOrds: see above.
%   G: updated G.
%
% Examples:
% coeff = graphfilt(G,s,g);
% coeff = graphfilt(G,s,g,cOrds);
% coeff = graphfilt(G,s,g,cOrds,opts_cOrds);
% coeff = graphfilt(G,s,g,[],opts_cOrds);
%
% Hamid Behjat

%-Checks. 
%--------------------------------------------------------------------------
if ~isfield(G,'N')
    G.N = size(G.A,1);
end
assert(size(s,1)==G.N,'G and G sigs mismatch in size.');
if ~isfield(G,'L') || isempty(G.L)
    G.L = sgwt_laplacian(G.A,'opt','normalized');
    G.lmax = sgwt_rough_lmax(G.L); % overwrite even if exist
end
if ~isfield(G,'lmax') || isempty(G.lmax)
    G.lmax = sgwt_rough_lmax(G.L);
end
arange = [0,G.lmax];

%-Estimate Chebyshev polynomial order to approximate kernels.
%--------------------------------------------------------------------------
if ~exist('opts_cOrds','var')
    opts_cOrds = [];
end
if ~exist('cOrds','var') || isempty(cOrds)
    [cOrds,~] = spgg_cheby_order_est(g,arange,opts_cOrds);
else
    assert(isstruct(cOrds),'cOrds should be a structure.');
    assert(isfield(cOrds,'kernel'));
    if ~opts_cOrds.justCheckKernels
        assert(isfield(cOrds,'tightframe'));
    end
end

%-Compute Chebyshev polynomial coefficients.
%--------------------------------------------------------------------------
c = cell(1,length(g));
for k = 1:length(g)
    if opts_cOrds.justCheckKernels
        d = cOrds.kernel(k);
    else
        d = cOrds.tightframe(k);
    end
    c{k} = sgwt_cheby_coeff(g{k},d,d+1,arange);
end

%-Decompose graph signals.
%--------------------------------------------------------------------------
coeff = sgwt_cheby_op(s,G.L,c,arange);
end