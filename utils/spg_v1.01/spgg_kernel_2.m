% added version '2' option for designing UMT filters such that the sum of 
% the abs values of the kernels forms a partition of unity rather than the 
% sum of the squared of the as values. 
% NOTE: maybe I should also change 'M' for this version so the uniformity
% constraint, cf. Eq.(28) in Behjat, et al. (2016) is satisfied. I have to
% look into that, but for the application I need this design now*, this is
% uniformity constraint is of no cencern. 
% * definig subdomains as benchmark for DILI.
% [25 June 2017]
%
% added (:) in:
% indNonEqual = [1; find(logical(diffOfEigs(:))) + 1];
% [13 June 2017]
%
% added option for binbased warping input
% [23 April 2017]
%
%
% added: 'fixForSparseSpectra' - for EMG project
% [17 Feb 2017]
%



%SPG_KERNEL_2 Second, third, ..., last spectral kernels of the uniform
% Meyer-type system of spectral kernels,
% as presented in (cf. Proposition 1):
%
% Behjat, et al., 'Signal-adapted tight frames on graphs',
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
%
% ------------------------------------------------------------
% Copyright (C) 2016, Hamid Behjat.
% This file is part of SPG (signal processing on graphs) package - v1.00
%
% Download: miplab.epfl.ch/software/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [r] = spgg_kernel_2(x,iScale,lmax,N_scales,varargin)

control_params={...
    'warping','none',...
    'E',[],...
    'subSampleWarping',0,...
    'fixForSparseSpectra',0,... % added in v1.01 [17 Feb 2017]
    'minSparseFactor',[],... % added in v1.01 [17 Feb 2017]
    'extrapolStep',[],...    % added in v1.01 [17 Feb 2017]
    'warpingBins',[],...     % added in v1.01 [23 April 2017]
    'vers',1}; % partition of unity formed by \sum|.|^2 [1], or \sum|.| [2] -- [25 June 2017]     

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if subSampleWarping==1 && fixForSparseSpectra==1
    error('This scenario has not been implemented.')
    % I think this even does not make sense to implement. [17 Feb 2017]
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

% NOTE 1
% This might appear as if its a bug since I am taking it to power of 2 for
% version 2 which is for forming partition of unity of abs values not
% squared abs values. But it is in fact correct.