function seqs = Ref2SeqSamples(Nbases,Nseq)
%% Tabulate homopolymers in human genome
%%
chrpath = 'C:\Users\MarkPratt\Documents\Genome\hg38\chroms';
chrStr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13', ...
    '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M'};
C = cell(size(chrStr));
N = zeros(size(chrStr));
parfor c = 1:numel(chrStr)
    fa = fastaread(fullfile(chrpath,['chr' chrStr{c} '.fa']));
    fa.Sequence = upper(fa.Sequence);
    fa.Sequence = fa.Sequence(fa.Sequence~='N');
    C{c}.Sequence = fa.Sequence;
    N(c) = numel(C{c}.Sequence);
end
S = zeros(1,sum(N),'uint8');
S = char(S);
ilast = 0;
for c = 1:numel(C)
    S(ilast+(1:N(c))) = C{c}.Sequence;
    C{c} = [];
    ilast = ilast + N(c);
end
    
seqs = zeros(Nbases,Nseq,'uint8');
seqs = char(seqs);
istart = randi(numel(S)-Nbases,Nseq,1);
ii = 0:Nbases-1;
for n = 1:Nseq
    seqs(:,n) = S(istart(n)+ii);
end
return