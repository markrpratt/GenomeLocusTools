function [R,ia,ib] = intersectGLS__(Ra,Rb)
%% function [R,ia,ib] = intersectGLS_(Ra,Rb)
% intersect two locus range sets confined to the same segment
% Ra & Rb are n x 2 arrays of start stop position on a single segment
% to do: return ia, ib indices for output, test coordinate intersection and
% bit logic v. current byte logic implementation
% warning: requires merged (and sorted) input
R  = [];
ia = [];
ib = [];

if isempty(Ra) || isempty(Rb)
    return;
elseif max(Ra(:,2)) < min(Rb(:,1))
    return;
elseif max(Rb(:,2)) < min(Ra(:,1))
    return;
end;

na = size(Ra,1);
nb = size(Rb,1);
rmax = max(na,nb);
R = zeros(rmax,2,'int32');
ia = zeros(rmax,1,'int32');
ib = zeros(rmax,1,'int32');
% Note: could be more intelligent about windowing to just plausible rows
ra = 1;
rb = 1;
r = 0;
while (1)
    % out of ranges
    if ra>na || rb>nb
        break;
    end;
    % disjoint
    if Ra(ra,2)<Rb(rb,1)
        ra = ra+1;
        continue;
    end;
    if Rb(rb,2)<Ra(ra,1)
        rb = rb+1;
        continue;
    end;
    % intersection
    r = r+1;
    R(r,1) = max(Ra(ra,1),Rb(rb,1));
    R(r,2) = min(Ra(ra,2),Rb(rb,2));
    ia(r) = ra;
    ib(r) = rb;
    switch sign(Ra(ra,2)-Rb(rb,2))      % who advances?
        case 1
            rb = rb + 1; 
        case -1
            ra = ra + 1;
        case 0
            ra = ra + 1;
            rb = rb + 1;
    end;
end;
R = R(1:r,:);
ia = ia(1:r);
ib = ib(1:r);
return;

