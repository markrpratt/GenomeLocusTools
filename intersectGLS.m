function [S,ia,ib] = intersectGLS(Sa,Sb,issorted)
%% function [S,ia,ib] = intersectGLS(Sa,Sb,issorted)
% intersect the locus set structures Sa and Sb returning merged S
% if issorted==true (default), operation should be faster and return ia,ib
% otherwise, a byte mask algo is used and ia, ib = []
% warning: default is to assume merged, sorted input
%%
if isempty(Sa) || isempty(Sb)
    S = []; ia = []; ib = [];
    return;
end;

if nargin < 3
    issorted = true;
end;

if isfield(Sa,'refID')
    S.refID = Sa.refID;
end;
%     if isfield(Sb,'refID')
%         if ~strcmp(Sa.refID,Sb.refID)
%             S = []; ia = []; ib = [];
%             return;
%         end;
%     end;
%     S.refID = Sa.refID;
% elseif isfield(Sb,'refID')
%     S.refID = Sb.refID;
% else
%     S.refID = [];
% end;

S.hdr = [Sa.hdr; Sb.hdr];
S.segNames = union(Sa.segNames, Sb.segNames);
[~,ia1,iSa] = intersect(S.segNames, Sa.segNames);
[~,ia2,iSb] = intersect(S.segNames, Sb.segNames);
I(ia1,1) = iSa;
I(ia2,2) = iSb;
ia = []; ib = [];
nR =0;
D(numel(S.segNames)).R = [];
D(end).ia = [];
D(end).ib = [];
parfor i = 1:numel(S.segNames)
    rra = find(Sa.R(:,1)==I(i,1));
    rrb = find(Sb.R(:,1)==I(i,2));
    if isempty(rra) || isempty(rrb)
        D(i).R = [];
        continue; 
    end;
    if issorted
        [R,D(i).ia,D(i).ib] = intersectGLS__(Sa.R(rra,2:3), Sb.R(rrb,2:3));
        D(i).ia = rra(D(i).ia);
        D(i).ib = rrb(D(i).ib);
    else
        R = intersectGLS_(Sa.R(rra,2:3), Sb.R(rrb,2:3));    % byte mask
    end;
    if isempty(R)
        D(i).R = [];
        D(i).ia = [];
        D(i).ib = [];
        continue;
    end;
    D(i).R(:,2:3) = R;
    D(i).R(:,1) = i;
    nR = nR + size(D(i).R,1);
end;

S.R = zeros(nR,3,'int32');
ia = zeros(nR,1,'int32');
ib = zeros(nR,1,'int32');
lastr = 0;
for i = 1:numel(S.segNames)
    nR = size(D(i).R,1);
    S.R(lastr+(1:nR),:) = D(i).R;
    ia(lastr+(1:nR)) = D(i).ia;
    ib(lastr+(1:nR)) = D(i).ib;
    lastr = lastr + nR;
end;
return;
