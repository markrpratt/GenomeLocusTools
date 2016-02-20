function S = mergeGLS(S1)
%% S = mergeGLS(S1)
% return compact and sorted locus set from input S1
% ranges sorted by start position and overlapping or abutting ranges are
% merged.
% to do: trace source indices in S1 of output s
%%
S.refID = S1.refID;
S.hdr = S1.hdr;
S.segNames =S1.segNames;
nR = 0;
for i = numel(S.segNames):-1:1
    D(i).R(:,2:3) = mergeGLS_(S1.R(S1.R(:,1)==i,2:3));
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
