function S = setdiffGLS(S1,S2)
%% function S = setdiffGLS(S1,S2)
% S = S1-S2, regions from S2 are removed from S1

S.refID = S1.refID;
S.hdr = [S1.hdr; S2.hdr];
S.segNames = union(S1.segNames, S2.segNames);
[~,ia1,is1] = intersect(S.segNames, S1.segNames);
[~,ia2,is2] = intersect(S.segNames, S2.segNames);
I(ia1,1) = is1;
I(ia2,2) = is2;
nR =0;
D.nR = 0;
D.R = [];
D(1:numel(S.segNames)) = D;
for i = 1:numel(S.segNames)
    D(i).nR = 0;
    if I(i,1)==0, continue; end;
    rr1 = find(S1.R(:,1)==I(i,1));
    rr2 = find(S2.R(:,1)==I(i,2));
    if isempty(rr2)
        D(i).R = S1.R(rr1,2:3);
        D(i).nR = numel(rr1);
        continue; 
    end;
    D(i).R = setdiffGLS_(S1.R(rr1,2:3), S2.R(rr2,2:3));
    D(i).nR = size(D(i).R,1);
end;
S.R = zeros(sum([D.nR]),3,'int32');
isSeg = [D.nR]>0;
S.segNames = S.segNames(isSeg);
D = D(isSeg);
lastr = 0;
for i = 1:numel(S.segNames)
    nR = size(D(i).R,1);
    S.R(lastr+(1:nR),2:3) = D(i).R;
    S.R(lastr+(1:nR),1) = i;
    lastr = lastr + nR;
end;
return;