function GLS = stdhg19ChrsGLS(inGLS)
%%
GLS.refID = inGLS.refID;
GLS.hdr = inGLS.hdr;
GLS.R = [];
GLS.segNames = {...
    'chr1'
    'chr2'
    'chr3'
    'chr4'
    'chr5'
    'chr6'
    'chr7'
    'chr8'
    'chr9'
    'chr10'
    'chr11'
    'chr12'
    'chr13'
    'chr14'
    'chr15'
    'chr16'
    'chr17'
    'chr18'
    'chr19'
    'chr20'
    'chr21'
    'chr22'
    'chrX'
    'chrY'
    'chrMT'};
if isempty(inGLS.R)
    return;
end;
[~,ii,jj] = intersect(GLS.segNames,inGLS.segNames);
I = zeros(size(inGLS.R(:,1)),'int32');
for n = 1:numel(jj)
    I(inGLS.R(:,1)==jj(n)) = ii(n);
end;
% I(I==0) = numel(GLS.segNames);
GLS.R = zeros(sum(I>0),3,'int32');
GLS.R(:,2:3) = inGLS.R(I>0,2:3);
GLS.R(:,1) = I(I>0);



