%% Tabulate homopolymers in human genome
clear;
chrpath = 'C:\Users\MarkPratt\Documents\Genome\hg38\chroms';
chrStr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13', ...
    '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M'};
hpll = 0:100;
hpl = zeros(numel(chrStr),4,numel(hpll));   % histogram of hpol length by nuc & chr
imap(int16('ACGTacgtNn')) = int8([1 2 3 4 1 2 3 4 0 0]);
tic;
parfor c = 1:numel(chrStr)
    fa = fastaread(fullfile(chrpath,['chr' chrStr{c} '.fa']));
    ifa = imap(int16(fa.Sequence));        % numeric sequence
    dn = find(diff(ifa)~=0);          % base transition positions
    lhp = diff(dn);              % hpol length between transitions
    ihp = ifa(dn(1:end-1)+1);
    for n = 1:4
        hpl(c,n,:) = hist(lhp(ihp==n),hpll);
    end
end
toc;
% Elapsed time is 181.168190 seconds.
