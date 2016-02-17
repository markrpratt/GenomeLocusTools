function R = setdiffGLS_(R1,R2)
%% function setdiffGLS_(R1,R2)
% subtract two locus range sets confined to the same segment
% R1 & R2 are n x 2 arrays of start stop position on a single segmentj
% Todo: implement with switch between list and logic vector algorithms
if isempty(R1)
    R  = [];
    return;
end;
if isempty(R2)
    R = R1;
    return;
end;

% determine range of starts and stops
minc = min([R1(:,1); R2(:,1)]);
maxc = max([R1(:,2); R2(:,2)]);
I = zeros(maxc-minc+1,1,'int8');

% set R1
for r = 1:size(R1,1)
    I((R1(r,1):R1(r,2))-minc+1) = 1;
end;
% unset R2
for r = 1:size(R2,1)
    I((R2(r,1):R2(r,2))-minc+1) = 0;
end;

if I(1) == 1
    iistart = int32([1; find(diff(I)>0)]);
else
    iistart = int32(find(diff(I)>0));
end;
if I(end) == 1
    iiend = int32([find(diff(I)<0); maxc]);
else
    iiend = int32(find(diff(I)<0));
end;
R = int32([iistart, iiend]+minc-1);
return;


