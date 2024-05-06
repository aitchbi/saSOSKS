function [E,G,g,cents,cOrds,opts_cOrds] = hb_get_coarse_eesd(G,S,Nu,sosksType,cOrds,opts_cOrds,varargin)
% HB_GET_COARSE_EESD computes a coarse representation of the ensemble
% energy spectral density (EESD) of a given graph signal set S defined on
% graph G. The computation is coarse, that is, at the resolution of
% slightly overlapping spectral bands as defined by a uniformly distributed
% system of spectral kernels (SOSKS); this input can instead be a desired
% system kernels given as a cell array. This is done by computing the
% ensemble amount of energy captured by each of the spectral kernels when
% they are applied to the signal set via Chebyshev polynomial
% approximation.

% Inputs:
%   G: graph strcuture, at minimum with field A (adjaceny matrix) 
%   S: graph signal set, one signal per column.
%   Nu: number of kernels in uniform system of kernels for eatimating graph
%   signal energies (default: 50). This input is ignored if sosksType is a
%   cell array specifying the sosks to be used.
%   sosksType: (optional) default: 'MeyerLike'. This input can also be a
%   cell array of function handles that specifies the SOSKS to be used.
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

d = inputParser;
addParameter(d,'ParsevalTolerance', 1e-2);
parse(d,varargin{:});
opts = d.Results;

%-Stuff.
%--------------------------------------------------------------------------
if ~exist('sosksType','var') || isempty(sosksType)
    sosksType = 'MeyerLike';
end
if iscell(sosksType)
    msg = 'cOrds should be input if sosksType is a cell array';
    assert(exist('cOrds','var'), msg);
    if exist('opts_cOrds','var')
        if isempty(opts_cOrds)
            opts_cOrds = struct;
            opts_cOrds.justCheckKernels = false;
        end
    else
        opts_cOrds = struct;
        opts_cOrds.justCheckKernels = false;
    end
else
    if ~exist('Nu','var') || isempty(Nu)
        Nu = 50;
    end
    if ~exist('cOrds','var') || isempty(cOrds)
        if exist('hb_get_precomp_cheby_ords.m','file')
            [cOrds,opts_cOrds] = hb_get_precomp_cheby_ords(Nu);
        else
            d1 = 'Function hb_get_precomp_cheby_ords.m not in path. ';
            d2 = 'Required Chebyshev orders will be computed.';
            warning([d1 d2]);
            cOrds = [];
            opts_cOrds = [];
        end
    end
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
if iscell(sosksType)
    cents = [];
    g = sosksType;
    Nu = length(g);
else
    assert(ischar(sosksType));
    switch sosksType
        case 'MeyerLike'
            sz = G.lmax/Nu;
            cents = (sz/2):(sz):(G.lmax-(sz/2));
            g = spgg_filter_design(G.lmax,Nu,'designtype','uniform_meyer_type');
            if 0
                % verify assigned kernel centers are reasonable
                figure; %#ok<UNRCH>
                hold on;
                for iK=1:Nu
                    e = 0:(G.lmax/(G.N-1)):G.lmax;
                    plot(e,g{iK}(e));
                    plot(cents(iK)*ones(1,2),[0 1],':k');
                end
            end

        case 'SplineType'
            error('extend code.')

        case 'sosks57'

        otherwise
            error('Unsupported kernel type for system of uniform kernels.')
    end
end

%-Compute coarse EESD.
%----------------------------------------------------------------------
E = zeros(1,Nu);
for k=1:N_gsigs
    s = S(:,k);
    s = s/norm(s);
    if isempty(cOrds)
        [coeffs,cOrds,G] = hb_graph_filt(G,s,g,cOrds,opts_cOrds);
    else
        [coeffs,~,G] = hb_graph_filt(G,s,g,cOrds,opts_cOrds);
    end
    E = E + cellfun(@norm,coeffs).^2;
    
end
E = E/N_gsigs;

assert(abs(sum(E)-1)<opts.ParsevalTolerance, ...
    'Fishy: the sum should be approximately 1.');
end
