function [xl,xh,pl,ph] = spgg_aux1_splinetype(so,sz,pou)
%
%
%
%
% Hamid Behjat
% Sep 2019.

[spl,~,~,x,ic,im] = hb_get_spline(so,sz,'approach','analytic');

if so==1
    xl = [];
    pl = [];
    xh = [];
    ph = [];
    
elseif ismember(so,[2,3])
    xl = cell(1,2);
    pl = cell(1,2);
    xh = cell(1,2);
    ph = cell(1,2);
    
    % first lower kernel
    xl{1} = x(ic:end);
    d1 = spl(ic:end);
    d2 = spl(im(2):end);
    pl{1} = d1+[d2;zeros(numel(d1)-numel(d2),1)];
    
    % 2nd lower kernel
    xl{2} = x(im(1):end);
    pl{2} = spl(im(1):end);
    
    % one to last upper kernel
    xh{1} = x(1:im(2));
    ph{1} = spl(1:im(2));
    
    % last upper kernel
    xh{2} = x(1:ic);
    d1 = spl(1:ic);
    d2 = spl(1:im(1));
    ph{2} = d1+[zeros(numel(d1)-numel(d2),1);d2];
    
    if strcmp(pou,'over2ndPower')
        pl{1} = sqrt(pl{1});
        pl{2} = sqrt(pl{2});
        ph{1} = sqrt(ph{1});
        ph{2} = sqrt(ph{2});
    end
    
elseif ismember(so,[4,5])
    xl = cell(1,3);
    pl = cell(1,3);
    xh = cell(1,3);
    ph = cell(1,3);
    
    % first lower kernel
    xl{1} = x(ic:end);
    d1 = spl(ic:end);
    d2 = spl(im(3):end);
    d3 = spl(im(4):end);
    pl{1} = d1...
        +[d2;zeros(numel(d1)-numel(d2),1)]...
        +[d3;zeros(numel(d1)-numel(d3),1)];
    
    % 2nd lower kernel
    xl{2} = x(im(2):end);
    pl{2} = spl(im(2):end);
    
    % 3rd lower kernel
    xl{3} = x(im(1):end);
    pl{3} = spl(im(1):end);
    
    % two to last upper kernel
    xh{1} = x(1:im(4));
    ph{1} = spl(1:im(4));
    
    % one to last upper kernel
    xh{2} = x(1:im(3));
    ph{2} = spl(1:im(3));
    
    % last upper kernel
    xh{3} = x(1:ic);
    d1 = spl(1:ic);
    d2 = spl(1:im(2));
    d3 = spl(1:im(1));
    ph{3} = d1...
        +[zeros(numel(d1)-numel(d2),1);d2]...
        +[zeros(numel(d1)-numel(d3),1);d3];
    
    if strcmp(pou,'over2ndPower')
        pl{1} = sqrt(pl{1});
        pl{2} = sqrt(pl{2});
        pl{3} = sqrt(pl{3});
        ph{1} = sqrt(ph{1});
        ph{2} = sqrt(ph{2});
        ph{3} = sqrt(ph{3});
    end
    
elseif ismember(so,[6,7])
    xl = cell(1,4);
    pl = cell(1,4);
    xh = cell(1,4);
    ph = cell(1,4);
    
    % first lower kernel
    xl{1} = x(ic:end);
    d1 = spl(ic:end);
    d2 = spl(im(4):end);
    d3 = spl(im(5):end);
    d4 = spl(im(6):end);
    pl{1} = d1...
        +[d2;zeros(numel(d1)-numel(d2),1)]...
        +[d3;zeros(numel(d1)-numel(d3),1)]...
        +[d4;zeros(numel(d1)-numel(d4),1)];
    
    % 2nd lower kernel
    xl{2} = x(im(3):end);
    pl{2} = spl(im(3):end);
    
    % 3rd lower kernel
    xl{3} = x(im(2):end);
    pl{3} = spl(im(2):end);
    
    % 4th lower end
    xl{4} = x(im(1):end);
    pl{4} = spl(im(1):end);
    
    % three to last kernel
    xh{1} = x(1:im(6));
    ph{1} = spl(1:im(6));
    
    % two to last upper kernel
    xh{2} = x(1:im(5));
    ph{2} = spl(1:im(5));
    
    % one to last upper kernel
    xh{3} = x(1:im(4));
    ph{3} = spl(1:im(4));
    
    % last upper kernel
    xh{4} = x(1:ic);
    d1 = spl(1:ic);
    d2 = spl(1:im(3));
    d3 = spl(1:im(2));
    d4 = spl(1:im(1));
    ph{4} = d1...
        +[zeros(numel(d1)-numel(d2),1);d2]...
        +[zeros(numel(d1)-numel(d3),1);d3]...
        +[zeros(numel(d1)-numel(d4),1);d4];
    
    if strcmp(pou,'over2ndPower')
        pl{1} = sqrt(pl{1});
        pl{2} = sqrt(pl{2});
        pl{3} = sqrt(pl{3});
        pl{4} = sqrt(pl{4});
        ph{1} = sqrt(ph{1});
        ph{2} = sqrt(ph{2});
        ph{3} = sqrt(ph{3});
        ph{4} = sqrt(ph{4});
    end
end