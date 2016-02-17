function R = intersectGLS_(R1,R2)
%% function R = intersectGLS_(R1,R2)
% intersect two locus range sets confined to the same segment
% R1 & R2 are n x 2 arrays of start stop position on a single segment
% to do: return ia, ib indices for output, test coordinate intersection and
% bit logic v. current byte logic implementation

if isempty(R1) || isempty(R2)
    R  = [];
    return;
end;
minc = min([R1(:,1); R2(:,1)]);
maxc = max([R1(:,2); R2(:,2)]);
I1 = zeros(maxc-minc+1,1,'int8');
I2 = zeros(maxc-minc+1,1,'int8');
% rstart = R1(r,1)-minc+1;
% rend = R1(r,2)-minc+1;
% parfor r = 1:size(R1,1)
%     I1(rstart(r):rend(r)) = 1;
% end;
for r = 1:size(R1,1)
    I1((R1(r,1):R1(r,2))-minc+1) = 1;
end;
for r = 1:size(R2,1)
    I2((R2(r,1):R2(r,2))-minc+1) = 1;
end;
I1 = I1 & I2;
if I1(1) == true
    iistart = int32([1; find(diff(I1)>0)+1]);
else
    iistart = int32(find(diff(I1)>0)+1);
end;
if I1(end) == 1
    iiend = int32([find(diff(I1)<0); maxc]);
else
    iiend = int32(find(diff(I1)<0));
end;
R = int32([iistart, iiend]+minc-1);
return;


