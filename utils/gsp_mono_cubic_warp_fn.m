% This function is an extended version of function gsp_mono_cubic_warp_fn.m
% from gspbox.sourceforge.net, which includes a number of additions and bug
% fixes as detailed in the changelog below:
%
% [28/06/2019] 
% Placed an initial check of the value 'cut'. With this
% check, I think the issue of [17/04/2017] gets resolved. As such, I
% commented out the [17/04/2017]; if in some scenraions the issue returns
% an error will be returned and I can then uncomment that fix.
%
% [29/01/2018] 
% added a little fix to prevent a bug, in step 2:
% 'Initialize tangents m at every data point'
%
% [29/01/2018] 
% added a sanity check at the end to make sure the 
% resulting interpolation is always monotic. For instance, for the 
% Alameda graph, interpolating at xxx = 0:0.001:G.lmax, leads to 
% nonmonotocin warp_ust at samples around 6400. 
%
% [24/01/2018] 
% added a check, line 79; the index sometines becomes 1.
%
% [17/04/2017]
% added a check, line 81; it sometimes goes out of bound. 


function [ interpolated_values ] = gsp_mono_cubic_warp_fn(x,y,x0)

cut=1e-4;
d = abs(min(diff(x0)));
if d<=cut %HB [28/06/2019]
    cut = d/10;
end

% Make sure data is sorted and monotonic
%
%   Url: http://gspbox.sourceforge.net/doc/filter_design/utils/gsp_mono_cubic_warp_fn.php

% Copyright (C) 2013 David I Shuman, Nathanael Perraudin.
% This file is part of GSPbox version 0.0.1
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
[x,x_ind]=sort(x,'ascend');
y=y(x_ind);
if ( isequal(sort(y,'ascend'),y)==0 && isequal(sort(y,'descend'),y)==0 )
    error('Data points are not monotonic');
end

% Monotonic cubic interpolation using the Fritsch-Carlson method
num_pts=length(x);
if length(y) ~= num_pts
    error('x and y vectors have different dimensions');
end

% 1. Compute slopes of secant lines
Delta=(y(2:end)-y(1:num_pts-1))./(x(2:end)-x(1:num_pts-1));

% 2. Initialize tangents m at every data point
m = (Delta(1:num_pts-2)+Delta(2:num_pts-1))/2;
m = m(:); % HB [29 Jan 2018]
m = [Delta(1);m;Delta(end)];

% 3. Check for equal y's to set slopes equal to zero
for k=1:num_pts-1
    if Delta(k)==0
        m(k)=0;
        m(k+1)=0;
    end
end

% 4. Initialize alpha and beta
alpha = m(1:num_pts-1)./Delta;
beta = m(2:num_pts)./Delta;

% 5. Make monotonic
for k=1:num_pts-1
    if alpha(k)^2+beta(k)^2 > 9
        tau=3/sqrt(alpha(k)^2+beta(k)^2);
        m(k)=tau*alpha(k)*Delta(k);
        m(k+1)=tau*beta(k)*Delta(k);
    end
end

% 6. Cubic interpolation
num_pts_to_interpolate=length(x0);
interpolated_values=zeros(size(x0));

for i=1:num_pts_to_interpolate
    [~,closest_ind]=min(abs(x-x0(i)));
    %if sign(x(closest_ind)-x0(i))<0 || ( sign(x(closest_ind)-x0(i))==0 && closest_ind < num_pts)
    if (x(closest_ind)-x0(i))<-cut || ( abs(x(closest_ind)-x0(i))<cut && closest_ind < num_pts)
        lower_ind=closest_ind;
    elseif closest_ind==1 % HB [24/01/2018]
        lower_ind = closest_ind;
    else
        lower_ind=closest_ind-1;
    end
    if lower_ind+1 > numel(x) % HB [17/04/2017]
        error('ooops.. replace back the code below, 3 lines..'); % HB [28/06/19]
        %fprintf('\n[gsp_mono_cubic.m] - HB: ooopps.. out of range.. skipping. \n')
        %interpolated_values(i) = NaN; % HB [28/06/19] Wouldn't it be better use the value next to it? like I do for when lower_ind==0
        %continue;
    %elseif lower_ind==0 % HB [28/06/19]
    %    lower_ind = 1;
    end
    h=x(lower_ind+1)-x(lower_ind);
    t=(x0(i)-x(lower_ind))/h;

    interpolated_values(i) = y(lower_ind)*(2*t^3-3*t^2+1) + h*m(lower_ind)*(t^3-2*t^2+t) + y(lower_ind+1)*(-2*t^3+3*t^2) + h*m(lower_ind+1)*(t^3-t^2); 
end

% sanity check to make sure result is monotonic & if not, to adjust it
if any(diff(interpolated_values)<0)
    d1 = interpolated_values;
    nd1 = numel(d1);
    d2 = find(diff(d1)<0,1);
    while ~isempty(d2)
        d3 = d2+2;
        chk = 1;
        while  chk && d1(d3)-d1(d2)<0
            if d3<nd1
                d3 = d3+1;
            else
                d3 = nd1;
                chk = 0;
            end
        end
        d4 = d3-d2;
        d5 = d1(d3)-d1(d2);
        d6 = d5/d4;
        d1(d2+1:d2+d4-1) = d1(d2)+d6:d6:d1(d3)-d6;
        d2 = find(diff(d1)<0,1);
    end
    if any(diff(d1)<0)
        error('fishy..');
    end
    interpolated_values = d1;
end
end


