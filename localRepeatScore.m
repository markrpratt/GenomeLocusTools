function [repScore, gcScore] = localRepeatScore(GLS,Nper, refpath)
%% function S = localRepeatScore(GLS,N)
% return repeat score for each region in GLS looking at periodicities 1:N
% S is matrix with N columns and one row for each entry in GLS
% each score in S is the average of forward and reverse shifted identity
% i.e. for period 2, S(i,2) = 0.5*(mean(ref(start:end)==ref((start:end)+2)
% + mean(ref(start:end)==ref((start:end)-2))
% gcscore is an additional vector of mean GC content since it was easy to
% compute here


if nargin < 3
    refpath = '/Users/mark/genome/hg19';
end;
segNames{24} = 'chrY';
segNames{23} = 'chrX';
for i = 1:22
    segNames{i} = sprintf('chr%d',i);
end;
repScore = zeros(size(GLS.R,1),Nper,'single');
gcScore = zeros(size(GLS.R,1),1,'single');

%%
for i = 1:numel(GLS.segNames)
    ss = find(GLS.R(:,1)==i);
    if isempty(ss), continue; end;
    S = faread(fullfile(refpath,['chr' segNames{i} '.fa']));
    Rseq = upper(S.Sequence);  
    clear S;
    R = GLS.R(ss,:);
    rs = zeros(numel(ss),Nper,'single');
    gc = zeros(numel(ss),1,'single');
    parfor j = 1:numel(ss)
        cc = R(j,2):R(j,3);
        for n = 1:N
            rs(j,n) = 0.5*( ...
                mean(Rseq(cc)==Rseq(cc+n)) + ...
                mean(Rseq(cc)==Rseq(cc-n)));
        end;
        gc(j) = mean(Rseq(cc)=='G' + Rseq(cc)=='C');
    end;
    repScore(ss,:) = rs;
    gcScore(ss) = gc;
end;
return;

