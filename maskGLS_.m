function Score = maskGLS_(inR,maskR)
%% function Score = maskGLS_(inR, maskR)
% return average  (fraction of overlap) of mask on inR records

if isempty(inR) || isempty(maskR)
    R  = [];
    return;
end;
minc = min([inR(:,1); maskR(:,1)]);
maxc = max([inR(:,2); maskR(:,2)]);
maskI = zeros(maxc-minc+1,1,'int8');

Score = zeros(size(inR,1),1,'single');
for r = 1:size(maskR,1)
    maskI((maskR(r,1):maskR(r,2))-minc+1) = 1;
end;
parfor r = 1:size(inR,1)
    Score(r) = mean(single(maskI((inR(r,1):inR(r,2))-minc+1)));
end;
return;
