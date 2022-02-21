function [w,e,w_ns,G] = hb_get_coarse_ee_warping(G,E,cents,L_sm)
% HB_GET_COARSE_EE_WARPING computes an energy equalizing (EE) warping that
% is required for building signal-adapted system of spectral kernels
% (saSOSKS). This warping is build based on the coarse estimate of the
% signal sets coarse EESD, obtained using hb_get_coare_eesd.m
%
% Inputs:
%   G: graph strcuture, at minimum with field A (adjaceny matrix) 
%   E: coarse EESD computed using hb_get_coare_eesd.m
%   cents: spectral center of kernels of the uniform SOSKS using for
%   obtaining E; output from hb_get_coare_eesd.m 
%   L_sm: (optional) length of moving avergage window for smoothing the
%   warping (default: 0.1*G.N, but it is recommended that you undersatnd
%   this parameter and set it yourself).
% 
% Outputs:
%   w   : warping 
%   e   : sudo eigs
%   w_ns: warping, non-smoothed
%   G    : updated G.
%
% Examples: 
% w = hb_get_coarse_ee_warping(G,E,cents);
% w = hb_get_coarse_ee_warping(G,E,cents,L_sm);
%
% Hamid Behjat

if ~exist('L_sm','var') || isempty(L_sm)
    L_sm = 0.1*G.N;
end
if ~isfield(G,'lmax') || isempty(G.lmax)
    G.lmax = sgwt_rough_lmax(sgwt_laplacian(G.A,'opt','normalized')); 
end

if sum(E)~=1
    E = E/sum(E);
end

C = cumsum(E);
C = [0 C];

e = 0:(G.lmax/(G.N-1)):G.lmax; % sudo eigs

if length(cents)==length(E)
    assert(cents(1)~=0,'fishy; center of 1st kernel is not 0.');
    cents = [0 cents];
else
     assert(length(cents)==(length(E)+1));
     assert(cents(1)==0);
end
cents(end) = G.lmax;

w_ns = interp1(cents,C,e);
w = smooth(w_ns,L_sm);

G.E = e(:);
end
