function S = BED2GLS(filename)
%% Read BED file into a Genome Locus Set structure
% translate n entries in BED file to a range matrix with rows
% [segID, start, stop] of int32
% segID are indices into S.segNames

S.refID = [];   % need to intelligently identify correct reference
S.hdr = [];
try 
    fdbed = fopen(filename);
catch fex
    disp(['unable to open ' filename]);
    return;
end;

nhdr = 0;
while 1
    hdrstr = fgets(fdbed);
%     if strcmp('chr',hdrstr(1:3)); % warning: chr not std on all refs
%         frewind(fdbed);
%         break;
%     end;
    if ~isempty(hdrstr) && hdrstr(1)~='#'  % start data - not comment or blank
        frewind(fdbed);
        break;
    end;
    nhdr = nhdr+1;
end;
    
S.hdr = cell(nhdr,1);
for i = 1:nhdr
    S.hdr{i} = fgets(fdbed);
end;

% minimum read
BEDdata = textscan(fdbed,'%s %d32 %d32 %*[^\n]', ...
    'delimiter','\t');
fclose(fdbed);
S.R = zeros(numel(BEDdata{1}),3,'int32');
[S.segNames, ~, S.R(:,1)] = unique(BEDdata{1});
S.R(:,2) = BEDdata{2}+1;
S.R(:,3) = BEDdata{3};
S.R = S.R(S.R(:,2)<=S.R(:,3),:);        % trim invalid ranges
% T = table(S.segNames(S.R(:,1)), S.R(:,2)-1, S.R(:,3), ...
%     'VariableNames',{'Seg','Start','End'});
% need to robustly extract additional bed fields
% full read 
% BEDdata = textscan(fdbed,'%s %d32 %d32 %s %f32 %c %*[^\n]', ...
%     'delimiter','\t');% 
% S.name = BEDdata{4};
% S.score = BEDdata{5};
% S.strand = BEDdata{6};

return;


