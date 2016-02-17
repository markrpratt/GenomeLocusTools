function L = dilateGLS(L,nbp,doMerge)
%% L = dilateGLS(L,nbp, doMerge)
% return locus set, L, with dilated ranges from input set L
% start and stop positions in each range will be extend by nbp
% default behavior is to merge returned locus set
% future: return indices of inputs mapped to output
% future: implement left, right and minimum size, other(?) options
if nargin < 3
    doMerge = true;
end;
L.R(:,2) = L.R(:,2)-nbp;
L.R(L.R(:,2)<1,2) = 1;
L.R(:,3) = L.R(:,3)+nbp;
if doMerge
    L = mergeGLS(L);
end;
return;
