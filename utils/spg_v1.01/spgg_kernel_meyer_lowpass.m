function y = spgg_kernel_meyer_lowpass(x,t,w,vers)
% SPGG_kernel_lowpass lowpass kernel using Meyer auxilary function, for a 
% given transition point and width w. A complementary highpass filter can
% be defiend (spgg_kernel_highpass.m) using the same parameters as used here
% in the lowpass filter. The resulting spectral kernel pair will have
% quadrature mirror property at the transition band.
%
% INPUTS:
%   x: position to return kernel value of. 
%   t: transition point
%   w: transition width. 
%   vers: 1 or 2. Partition of unity type. 
%        1: \sum|.|^2 
%        2: \sum|.| 
%        vers=2 can be used, for instance, in definig subdomains as 
%        benchmark for DIBSI. 
%
% OUTPUTS:
%   y: kernel values at points x. 
%
% Based on: spgg_kernel_1.m
%
% Also see:
% spgg_kernel_highpass.m
%
%
% Hamid Behjat
% Nov 2020.

if ~exist('vers','var') || isempty(vers)
    vers = 1; 
end

M = 2*w/(2*t-w)+1; 
q = 1/(M-1);
a = t-w/2;

y = zeros(size(x));

y(x>=0 & x<=a) = 1;

d = find(x>a & x<=M*a);
y(d) = cos(pi/2*(meyeraux(q*(x(d)/a -1))));

if vers==2
    y = y.^2; % see NOTE. 
end
end

function y = meyeraux(x)
p = [-20 70 -84 35 0 0 0 0];
y = polyval(p,x);
end

% NOTE
% This might appear as if its a bug since I am taking it to power of 2 for 
% version 2 which is for forming partition of unity of abs values not 
% squared abs values. But it is correct. 


