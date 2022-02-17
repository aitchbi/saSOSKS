function [d3,d4] = spgg_aux2_splinetype(so,sz,iN,N,sc,lmax,pou,xl,xh,pl,ph)
%
%
%
%
% Hamid Behjat
% Sep 2019.

[spl,~,~,x] = hb_get_spline(so,sz,'approach','analytic');
splsqrt = sqrt(spl);


if so==1
    d1 = sc*(x+iN-1);
    d2 = splsqrt;
    
elseif ismember(so,[2,3])
    switch iN
        case 1
            d1 = sc*(xl{1}+iN-1);
            d2 = pl{1};
        case 2
            d1 = sc*(xl{2}+iN-1);
            d2 = pl{2};
        case N-1
            d1 = sc*(xh{1}+iN-1);
            d2 = ph{1};
        case N
            d1 = sc*(xh{2}+iN-1);
            d2 = ph{2};
        otherwise
            d1 = sc*(x+iN-1);
            if strcmp(pou,'over2ndPower')
                d2 = splsqrt;
            else
                d2 = spl;
            end
    end
    
elseif ismember(so,[4,5])
    switch iN
        case 1
            d1 = sc*(xl{1}+iN-1);
            d2 = pl{1};
        case 2
            d1 = sc*(xl{2}+iN-1);
            d2 = pl{2};
        case 3
            d1 = sc*(xl{3}+iN-1);
            d2 = pl{3};
        case N-2
            d1 = sc*(xh{1}+iN-1);
            d2 = ph{1};
        case N-1
            d1 = sc*(xh{2}+iN-1);
            d2 = ph{2};
        case N
            d1 = sc*(xh{3}+iN-1);
            d2 = ph{3};
        otherwise
            d1 = sc*(x+iN-1);
            if strcmp(pou,'over2ndPower')
                d2 = splsqrt;
            else
                d2 = spl;
            end
    end
    
elseif ismember(so,[6,7])
    switch iN
        case 1
            d1 = sc*(xl{1}+iN-1);
            d2 = pl{1};
        case 2
            d1 = sc*(xl{2}+iN-1);
            d2 = pl{2};
        case 3
            d1 = sc*(xl{3}+iN-1);
            d2 = pl{3};
        case 4
            d1 = sc*(xl{4}+iN-1);
            d2 = pl{4};
        case N-3
            d1 = sc*(xh{1}+iN-1);
            d2 = ph{1};
        case N-2
            d1 = sc*(xh{2}+iN-1);
            d2 = ph{2};
        case N-1
            d1 = sc*(xh{3}+iN-1);
            d2 = ph{3};
        case N
            d1 = sc*(xh{4}+iN-1);
            d2 = ph{4};
        otherwise
            d1 = sc*(x+iN-1);
            if strcmp(pou,'over2ndPower')
                d2 = splsqrt;
            else
                d2 = spl;
            end
    end
end

[d3,d4] = hb_pvt(d1,d2,lmax);

if length(unique(d3(:)))~=length(d3) % 27 June 2019 --- d3 must be unique, othewise interp1 will run into error when calling g{.}
    %error('something is fishy..')
    d5 = [true;logical(diff(d3))];
    d3 = d3(d5);
    d4 = d4(d5);
end
end

function [d3,d4] = hb_pvt(d1,d2,lmax)
d1 = d1(:);
d2 = d2(:);
if d1(1)==0
    d3 = [d1(:);lmax];
    d4 = [d2;0];
elseif d1(end)==lmax
    d3 = [0;d1];
    d4 = [0;d2];
else
    d3 = [0;d1;lmax];
    d4 = [0;d2;0];
end
end