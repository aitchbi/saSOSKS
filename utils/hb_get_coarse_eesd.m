function [E,G,g_unif,cents,cOrds,opts_cOrds] = hb_get_coarse_eesd(G,S,Nu,unifType,cOrds,opts_cOrds)
% HB_GET_COARSE_EESD computes a coarse representation of the ensemble
% energy spectral density (EESD) of a given graph signal set S defined on
% graph G. The computation is coarse, that is, at the resolution of
% slightly overlapping spectral bands as defined by a uniformly distributed
% system of spectral kernels. This is done by computing the ensemble amount
% of energy captured by each of the spectral kernels when they are applied
% to the signal set via Chebyshev polynomial approximation.

% Inputs:
%   G: graph strcuture, at minimum with field A (adjaceny matrix) 
%   S: graph signal set, one signal per column.
%   Nu: number of kernels in uniform system of kernels for eatimating graph
%   signal energies (default: 50).
%   unifType: (optional) default: 'MeyerLike'.
%   cOrds: (optional) Chebyshev polynomial orders for applying the uniform
%   system of spectral kernels to the signal to estimate spectral energy
%   density.
%   opts_cOrds: (optional) polynomial order estimation parameters, min/max
%   order, etc. 
%
% Outputs:
%   E: estimated spectral energy density.
%   G    : updated G.
%   g_unif: uniform SOSKS ued for estimating spectral energy density.
%   cents: spectral centers of uniform SOSKS.
% 
% Examples: 
% E = hb_get_coarse_eesd(G,S,Nu,[],cOrds,opts_cOrds);
%
% Hamid Behjat

%-Stuff.
%--------------------------------------------------------------------------
if ~exist('Nu','var') || isempty(Nu)
    Nu = 50;
end
if ~exist('unifType','var') || isempty(unifType)
    unifType = 'MeyerLike';
end
if ~exist('cOrds','var') || isempty(cOrds)
    if exist('hb_get_precomp_cheby_ords.m','file')
        [cOrds,opts_cOrds] = hb_get_precomp_cheby_ords(Nu);
    else
        d1 = 'Function hb_get_precomp_cheby_ords.m not in path. ';
        d2 = 'Required Chebyshev orders will be computed.';
        warning([d1 d2]);
        cOrds = [];
    end
end
if ~exist('opts_cOrds','var')
    opts_cOrds = [];
end
N_gsigs = size(S,2);
if ~isfield(G,'L') || isempty(G.L)
    G.L = sgwt_laplacian(G.A,'opt','normalized');
    G.lmax = sgwt_rough_lmax(G.L); % overwrite even if exist
end
if ~isfield(G,'lmax') || isempty(G.lmax)
    G.lmax = sgwt_rough_lmax(G.L);
end

%-Get uniform SOSKS
%--------------------------------------------------------------------------
switch unifType
    case 'MeyerLike'
        sz = G.lmax/Nu;
        cents = (sz/2):(sz):(G.lmax-(sz/2));
        g_unif = spgg_filter_design(G.lmax,Nu,'designtype','uniform_meyer_type');
        if 0
            % verify assigned kernel centers are reasonable
            figure; %#ok<UNRCH>
            hold on;
            for iK=1:Nu
                e = 0:(G.lmax/(G.N-1)):G.lmax;
                plot(e,g_unif{iK}(e));
                plot(cents(iK)*ones(1,2),[0 1],':k');
            end
        end
    case 'SplineType'
        error('extend code.')
    otherwise
        error('Unsupported kernel type for system of uniform kernels.')
end

%-Compute coarse EESD.
%----------------------------------------------------------------------
E = zeros(1,Nu);
for k=1:N_gsigs
    s = S(:,k);
    s = s/norm(s);
    if isempty(cOrds)
        [coeffs,cOrds,G] = hb_graph_filt(G,s,g_unif,cOrds,opts_cOrds);
    else
        [coeffs,~,G] = hb_graph_filt(G,s,g_unif,cOrds,opts_cOrds);
    end
    E = E + cellfun(@norm,coeffs).^2;
    
end
E = E/N_gsigs;
end
