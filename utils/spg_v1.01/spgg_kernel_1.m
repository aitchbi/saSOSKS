function [r] = spgg_kernel_1(x,lmax,N_scales,varargin)

control_params={...
    'warping','none',...
    'E',[],...
    'subSampleWarping',0,...
    'fixForSparseSpectra',0,... 
    'minSparseFactor',[],... 
    'extrapolStep',[],...    
    'warpingBins',[],...     
    'vers',1}; % partition of unity formed by \sum|.|^2 [1], or \sum|.| [2]

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if subSampleWarping==1 && fixForSparseSpectra==1
    error('This scenario has not been implemented.')
    % HB: this even may not make sense to implement.
end

% Apply the Warping to x (i.e. x = w(x))
x = spgg_apply_warping(x,warping,lmax,E,warpingBins,subSampleWarping,...
    fixForSparseSpectra,minSparseFactor,extrapolStep);

t = 1;
M = 2.73; 
q = 1/(M-1);
a = lmax/(M*N_scales-N_scales+2);
r = zeros(size(x));

r1 = find(x>=0 & x<=a);
dummy = ones(size(x));
r(r1) = dummy(r1);

r2 = find(x>a & x<=M*a);
r(r2) = cos(pi/2*(meyeraux(q*(t*(x(r2))/a -1))));

if vers==2
    r = r.^2; % NOTE 1 
end

end

function y = meyeraux(x)
p = [-20 70 -84 35 0 0 0 0];
y = polyval(p,x);
end

% NOTE 1 This might appear as if its a bug since it's taken to power 2 for
% version 2, which is for forming partition of unity of abs values not
% squared abs values. But it is correct.


