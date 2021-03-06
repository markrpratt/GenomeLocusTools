function S = BED2GLS(filename)
%% Read BED file into a Genome Locus Set structure
% translate n entries in BED file to a range matrix with rows
% [segID, start, stop] of int32
% segID are indices into S.segNames
% warning: segment names are stripped of any preceding 'chr'
% to do: read additional fields on deman

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
% strip chr prefix
S.segNames = regexprep(S.segNames,'chr','');
S.segNames = regexprep(S.segNames,'MT','M');
% ischrseg = strncmp('chr',S.segNames,3);
% for n = 1:numel(S.segNames)
%     if ischrseg(n)
%         S.segNames{n} = S.segNames{n}(4:end);
%     end;
% end;

StdSegNames = { ...
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', ...
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M'};
[~,ii,jj] = intersect(StdSegNames,S.segNames);
kk = setdiff(1:numel(S.segNames),jj);  % move non-std names in GLS to end
I = zeros(size(S.R,1),1,'int32');
for n = numel(kk):-1:1
   NewSegNames{n+numel(StdSegNames)} = S.segNames{kk(n)};
    I(S.R(:,1)==kk(n)) = n+numel(StdSegNames);
end;
% populate first 25 indices with standard chrs
NewSegNames(1:numel(StdSegNames)) = StdSegNames;
for n = 1:numel(jj)  
    I(S.R(:,1)==jj(n),1) = ii(n);
end;
S.segNames = NewSegNames;
S.R(:,1) = I;
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


