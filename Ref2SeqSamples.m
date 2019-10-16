function [jseq,useq] = Ref2SeqSamples(Nbases,Nseq)
%% produce table of sampled seqences from human reference
% jseq - base array (Nbases x Nseq) in uint8 indices
% useq - base symbol array
%%
chrpath = 'C:\Users\MarkPratt\Documents\Genome\hg38\chroms';
chrStr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13', ...
    '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M'};
C = cell(size(chrStr));
N = zeros(size(chrStr));
useq = {'A','C','G','T'};
parfor c = 1:numel(chrStr)
    fa = fastaread(fullfile(chrpath,['chr' chrStr{c} '.fa']));
    fa.Sequence = upper(fa.Sequence);
    fa.Sequence = fa.Sequence(fa.Sequence~='N');
    C{c}.jseq = uint8(size(fa.Sequence));
    for i = 1:numel(useq)
        C{c}.jseq(fa.Sequence==useq{i}) = uint8(i);
    end
    N(c) = numel(C{c}.jseq);
end

% create contiguous sequence
S = zeros(1,sum(N),'uint8');
ilast = 0;
for c = 1:numel(C)
    S(ilast+(1:N(c))) = C{c}.jseq;
    C{c} = [];
    ilast = ilast + N(c);
end
clear C; 
jseq = zeros(Nbases,Nseq,'uint8');
istart = randi(numel(S)-Nbases,Nseq,1);
ii = 0:Nbases-1;
parfor n = 1:Nseq
    jseq(:,n) = S(istart(n)+ii);
end

return