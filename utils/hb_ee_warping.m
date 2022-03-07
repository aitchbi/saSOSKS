function [w,e,w_ns,G,g_unif,E,cents] = hb_ee_warping(G,S,L_sm,Nu,unifType,cOrds,opts_cOrds)
% HB_EE_WARPING computes an energy equalizing warping required for building
% signal-adapted system of spectral kernels based on the energy content of
% the given sample signal set S.
%
% Inputs:
%   G: graph strcuture, at minimum with field A (adjaceny matrix) 
%   S: graph signal set, one signal per column.
%   cOrds: (optional) Chebyshev polynomial orders for applying the uniform
%   system of spectral kernels to the signal to estimate energy content.
%   L_sm: length of moving avergage window for smoothing the warping
%   (default: 0.1*G.N).
%   Nu: number of kernels in uniform system of kernels for eatimating graph
%   signal energies (default: 50).
% 
% Outputs:
%   w   : warping 
%   e   : sudo eigs
%   w_ns: warping, non-smoothed
%   G    : updated G.
%
% Examples: 
% w = hb_ee_warping(G,S);
% w = hb_ee_warping(G,S,L_sm);
% w = hb_ee_warping(G,S,L_sm,Nu);
% w = hb_ee_warping(G,S,L_sm,Nu);
% w = hb_ee_warping(G,S,L_sm,Nu,unifType);
% [w,e,~,G] = hb_ee_warping(G,S,[],Nu);
%
% Hamid Behjat

%-Stuff.
%--------------------------------------------------------------------------
if ~exist('L_sm','var') || isempty(L_sm)
    L_sm = 0.1*G.N;
end
if ~exist('Nu','var') || isempty(Nu)
    Nu = 50;
end
if ~exist('unifType','var') || isempty(unifType)
    unifType = 'MeyerLike';
end
if ~exist('cOrds','var')
    cOrds = [];
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

%-Get uniform Meyer-like SOSKS.
%--------------------------------------------------------------------------
g_unif = spgg_filter_design(G.lmax,Nu,'designtype','uniform_meyer_type');

%-Estimate signal set's ensemble energy.
%--------------------------------------------------------------------------
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
C = cumsum(E);
C = [0 C];

%-Approximate center of Uniform Meyer-like kernels.
%--------------------------------------------------------------------------
switch unifType
    case 'MeyerLike'
        sz = G.lmax/Nu;
        cents = (sz/2):(sz):(G.lmax-(sz/2));
        cents = [0 cents];
        cents(end) = G.lmax;
        if 0
            % verify assigned kernel centers are reasonable
            figure; %#ok<UNRCH>
            hold on;
            for iK=1:Nu
                plot(e,g_unif{iK}(e));
                plot(cents(iK+1)*ones(1,2),[0 1],':k');
            end
        end
    case 'SplineType'
        error('extend code.')
    otherwise
        error('Unsupported kernel type for system of uniform kernels.')
end

%-Upsample C_energy to get warping.
%--------------------------------------------------------------------------
e = 0:(G.lmax/(G.N-1)):G.lmax; % sudo eigs
w_ns = interp1(cents,C,e);
w = smooth(w_ns,L_sm);

% update G.
G.E = e(:);
cents = cents(2:end);
end
