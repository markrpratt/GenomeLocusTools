function R = mergeGLS_(inR)
%% function R = mergeGLS_(R)
% merge a locus set confined to a single segment
% R is an n x 2 array of start stop position on a single segment
R = sortrows(inR);
openr = 1;                                  % index of open interval
for r = 2:size(R,1)   
    if R(r,1)<=R(openr,2)+1   % overlap or bookend - merge with open interval
        R(openr,2) = max(R([openr; r],2));  % expand open range if nec.
        R(r,:) = 0;                         % clear current row
    else
        openr = r;
    end;
end;
R = R(R(:,1)>0,:);
return; 


