function Score = maskGLS(inS,maskS)
%% function Score = maskGLS(S,Smask)
% intersect the locus set structures S and Smask returning the average of
% mask on every record of S

%%
if isempty(inS)
    Score = []; 
    return;
end;
Score = zeros(size(inS.R,1),1,'single');
if isempty(maskS)
    return;
end;

segNames = union(inS.segNames, maskS.segNames);
[~,ia1,iSa] = intersect(segNames, inS.segNames);
[~,ia2,iSb] = intersect(segNames, maskS.segNames);
I(ia1,1) = iSa;
I(ia2,2) = iSb;
for i = 1:numel(segNames)
    rrin = find(inS.R(:,1)==I(i,1));
    rrmask = find(maskS.R(:,1)==I(i,2));
    if isempty(rrin) || isempty(rrmask), continue; end;
    Score(rrin) = maskGLS_(inS.R(rrin,2:3), ...
        maskS.R(rrmask,2:3));    % byte mask
end;

return;
