function GLS = GRCtoStdChrSegNames(inGLS)
%%
segNameLookup = {...
    '1', 'chr1'
    '2'  'chr2'
    '3'  'chr3'
    '4'  'chr4'
    '5'  'chr5'
    '6', 'chr6'
    '7', 'chr7'
    '8', 'chr8'
    '9', 'chr9'
    '10', 'chr10'
    '11', 'chr11'
    '12', 'chr12'
    '13', 'chr13'
    '14', 'chr14'
    '15', 'chr15'
    '16', 'chr16'
    '17', 'chr17'
    '18', 'chr18'
    '19', 'chr19'
    '20', 'chr20'
    '21', 'chr21'
    '22', 'chr22'
    'X', 'chrX'
    'Y', 'chrY'
    'MT', 'chrMT'};
GLS.refID = inGLS.refID;
GLS.hdr = inGLS.hdr;
GLS.R = [];
GLS.segNames = segNameLookup(:,2);
if isempty(inGLS.R)
    return;
end;
[~,ii,jj] = intersect(segNameLookup(:,1),inGLS.segNames);
I = zeros(size(inGLS.R(:,1)),'int32');
for n = 1:numel(jj)
    I(inGLS.R(:,1)==jj(n)) = ii(n);
end;
% I(I==0) = numel(GLS.segNames);
GLS.R = zeros(sum(I>0),3,'int32');
GLS.R(:,2:3) = inGLS.R(I>0,2:3);
GLS.R(:,1) = I(I>0);



