%% Test GenomeLocusSet tools

clear;
nlines = 20;
S0.refID = 'hg19';
S0.hdr = [];
S0.segNames = { ...
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', ...
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT'};
S0.R(:,2:3) = sort(randi(1e6,nlines,2),2);
S0.R(:,1) = randi(25,nlines,1);
GLS2BED(S0,'test0.bed');
S1 = mergeGLS(S0);
GLS2BED(S1,'test1.bed');
S2 = intersectGLS(S1,S1);
S3 = setdiffGLS(S2,S1);
statsGLS(S3,true);