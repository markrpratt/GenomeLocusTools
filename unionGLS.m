function S = unionGLS(S1,S2,merge)
if nargin < 3
    merge = true;
end;
    
%%
S.refID = S1.refID;
S.hdr = [S1.hdr; S2.hdr];
S.segNames = union(S1.segNames, S2.segNames);
[~,ia1,is1] = intersect(S.segNames, S1.segNames);
[~,ia2,is2] = intersect(S.segNames, S2.segNames);
I(ia1,1) = int32(is1);
I(ia2,2) = int32(is2);
nR =0;

for i = numel(S.segNames):-1:1
    if merge
        D(i).R(:,2:3) = mergeGLS_([ ...
            S1.R(S1.R(:,1)==I(i,1),2:3);
            S2.R(S2.R(:,1)==I(i,2),2:3)]);
    else
        D(i).R(:,2:3) = [ ...
            S1.R(S1.R(:,1)==I(i,1),2:3);
            S2.R(S2.R(:,1)==I(i,2),2:3)];
    end;
    D(i).R(:,1) = i;
    nR = nR + size(D(i).R,1);
end;

S.R = zeros(nR,3,'int32');
lastr = 0;
for i = 1:numel(S.segNames)
    nR = size(D(i).R,1);
    S.R(lastr+(1:nR),:) = D(i).R;
    lastr = lastr + nR;
end;
