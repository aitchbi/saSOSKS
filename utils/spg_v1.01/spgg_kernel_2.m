function [r] = spgg_kernel_2(x,iScale,lmax,N_scales,varargin)

control_params={...
    'warping','none',...
    'E',[],...
    'subSampleWarping',0,...
    'fixForSparseSpectra',0,... 
    'minSparseFactor',[],... 
    'extrapolStep',[],...    
    'warpingBins',[],...     
    'vers',1}; % partition of unity formed by \sum|.|^2 [1], or \sum|.| [2] -- [25 June 2017]     

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if subSampleWarping==1 && fixForSparseSpectra==1
    error('This scenario has not been implemented.')
    % HB: this may even not make sense to implement.
end

if isempty(minSparseFactor) %#ok<NODEF>
    minSparseFactor = lmax/10;
    if isempty(extrapolStep) %#ok<NODEF>
        extrapolStep = minSparseFactor/2;
    end
end

t = 1;
M = 2.73;
q= 1/(M-1);
a = lmax/(M*N_scales-N_scales+2);

r=zeros(size(x));

Delta = M*a-a;
%scaleShift = Delta*(iScale-1);

% Apply the Warping to x (i.e. x = w(x))
x = spgg_apply_warping(x,warping,lmax,E,warpingBins,...
    subSampleWarping,fixForSparseSpectra,minSparseFactor,extrapolStep);

r1 = find(x>= a+Delta*(iScale-1) & x<= M*a + Delta*(iScale-1));
r(r1) = sin(pi/2*(meyeraux(q*(t*(x(r1)- Delta*(iScale-1))/a -1))));

if iScale<N_scales
    r2 = find(x> M*a+Delta*(iScale-1) & x <= M*a + Delta*iScale);
    r(r2) = cos(pi/2*(meyeraux(q*(t*(x(r2)- Delta*iScale)/a -1))));
else
    r2 = find(x> M*a+Delta*(iScale-1));
    dummy = ones(size(x));
    r(r2) = dummy(r2);
end

if vers==2
    r = r.^2; % NOTE 1
end
end

function y = meyeraux(x)
p = [-20 70 -84 35 0 0 0 0];
y = polyval(p,x);
end

% NOTE 1 This might appear as if its a bug since I am taking it to power of
% 2 for version 2, which is for forming partition of unity of abs values
% not squared abs values. But it is correct.