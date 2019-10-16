function GLS2BED(S,filename,mode,name,score)
%% Read BED file into a Genome Locus Set structure
% translate n entries in BED file to a range matrix with rows
% [segID, start, stop] of int32
% segID are indices into S.segNames
if nargin < 3
    mode = 'w';
end

try 
    fdbed = fopen(filename,mode);
catch fex
    disp(['unable to open ' filename]);
    return;
end

if strcmp(mode,'w')
    for i = 1:numel(S.hdr)
        fprintf(fdbed,'%s\n',S.hdr{i});
    end
end

if ~isempty(S.refID)
    fprintf(fdbed, '# refID: %s\n', S.refID);
end
for i = 1:size(S.R,1)
    fprintf(fdbed,'%s\t%d\t%d\n', S.segNames{S.R(i,1)}, S.R(i,2)-1, S.R(i,3));
end
fclose(fdbed);

return;


