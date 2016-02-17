function L = erodeGLS(L,nbp,minsize)
%% L = erodeGLS(L,nbp, minsize)
% return locus set, L, with eroded ranges from input set L
% start and stop positions in each range will be reduced by nbp
% if the resulting range is less than minsize (default=0), the range will
% be expanded to be at least minsize

if nargin < 4
    minsize = 0;
end;
L.R(:,2) = L.R(:,2)+nbp;
L.R(:,3) = L.R(:,3)-nbp;

% remove negative and zero size ranges
if minsize == 0
    L.R = L.R(L.R(:,2)<=L.R(:,3),:);
    return;
end;

l = diff(L.R(:,2:3),1,2)+1;
d = single(max(0,minsize-l));
% expand to minsize
L.R(:,2) = L.R(:,2)-int32(floor(d/2));
L.R(:,3) = L.R(:,3)+int32(ceil(d/2));
L = mergeGLS(L);

return;