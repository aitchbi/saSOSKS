function [w, e, w_ns, Y] = hb_get_consensus_ee_warping(A,B,varargin)
% HB_GET_CONSENSUS_EE_WARPING builds a consensus energy-equalizing warping
% that can be used to build signal-adapted system of spectral graph kernels
% (saSOSKS). The consensus related to the that warping is defined based on
% varying graphs (and graphs signals) that do not need to be of the same
% size. 
%
% Inputs:
%   A: 
%   B: 
%
% Outputs: 
%
%
% Dependencies: 
% 
% Hamid Behajt 

d = inputParser;
addParameter(d,'SmoothSpan', []);
addParameter(d,'TypeB', 'energy'); % 'GFT' | 'energy'
addParameter(d,'InterpolationStep', []);
addParameter(d,'NumberOfInterpolationPoints', []);
parse(d,varargin{:});
opts = d.Results;

chkopts(opts);
    
[A1, A2, B1, B2] = verifyab(A,B);

K = size(A1,1);

if isempty(A2)
    TwoSets = 0;
else
    TwoSets = 1;
end

for k=1:K
    A1{k}    = abs(A1{k}); % all eigs positive
    A1{k}(1) = 0;         % first eig 0
    if TwoSets
        A2{k}    = abs(A2{k});
        A2{k}(1) = 0;
    end
end

switch opts.TypeB
    case 'GFT'
        % convert GFT coeffs to energy
        B1 = cellfun(@(x) x.^2, B1, 'UniformOutput', false); 
        if TwoSets
            B2 = cellfun(@(x) x.^2, B2, 'UniformOutput', false);
        end
    case 'energy'
        % B already is energy 
end

lmax = min(cellfun(@max, [A1; A2]));

if isempty(opts.InterpolationStep)
    if isempty(opts.NumberOfInterpolationPoints)
        d = cellfun(@length, [A1; A2]);
        M = round(mean(d(:)));
    else
        M = opts.NumberOfInterpolationPoints;
    end
else
    M = round(lmax/opts.InterpolationStep);
end
e = linspace(0,lmax,M);

BC1 = cellfun(@(x) cumsum(x,1), B1, 'UniformOutput', false);

if TwoSets
    BC2 = cellfun(@(x) cumsum(x,1), B2, 'UniformOutput', false);
end

Y_sum = zeros(M,1);
Y     = cell(K,1);
for k=1:K
    J = size(B1{k},2);
    Y{k} = zeros(M,J);
    A1k = A1{k}(:);
    if TwoSets
        A2k = A2{k}(:);
    end
    for j=1:J
        d = interp1(A1k, BC1{k}(:,j), e);
        if TwoSets
            d = abs(d - interp1(A2k, BC2{k}(:,j), e));
            d = cumsum(d);
        end
        Y{k}(:,j) = d;
    end
    if k==1
        J_sum = J;
        Y_sum = sum(Y{k},2);
    else
        J_sum = J_sum + J;
        Y_sum = Y_sum + sum(Y{k},2);
    end
end
w_ns = Y_sum/J_sum;
w_ns = w_ns./w_ns(end);

%-Smooth.
SS = opts.SmoothSpan;
if isempty(SS)
    SS = 0.02*lmax;
else
    assert(SS<1);
end
L = (SS/lmax)*M;
w = smooth(w_ns, L);
w = w./w(end);
end

%==========================================================================
function [A1, A2, B1, B2] = verifyab(A,B)
% A: Nx1 or Nx2 cell array. 
% B: Nx1 or Nx2 cell array.
% Both should be either a vector or matrix.
assert(or(isvector(A), ismatrix(A)));
assert(or(isvector(B), ismatrix(B)));
if isvector(A)
    assert(isvector(B));
    TwoSets = 0;
else
    TwoSets = 1;
end
N = size(A,1);
assert(size(B,1)==N);
for k=1:N
    assert(size(A{k,1},1)==size(B{k,1},1)); % lh
    if TwoSets
        assert(size(A{k,2},1)==size(B{k,2},1)); % rh
    end
end
if TwoSets
    A1 = A(:,1);
    B1 = B(:,1);
    A2 = A(:,2);
    B2 = B(:,2);
else
    A1 = A;
    B1 = B;
    A2 = [];
    B2 = [];
end
end

%==========================================================================
function chkopts(opts)
% at least one of them should be empty
if isempty(opts.InterpolationStep)
    return;
else
    assert(isempty(opts.NumberOfInterpolationPoints));
end
end