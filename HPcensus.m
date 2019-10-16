%% Tabulate homopolymers in human genome
clear;
chrpath = 'C:\Users\MarkPratt\Documents\Genome\hg38\chroms';
chrStr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13', ...
    '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M'};
hpll = 0:100;
hpl = zeros(numel(chrStr),4,numel(hpll));   % histogram of hpol length by nuc & chr
imap(int16('ACGTacgtNn')) = int8([1 2 3 4 1 2 3 4 0 0]);
tic;
Lchr = zeros(size(chrStr));
parfor c = 1:numel(chrStr)
    fa = fastaread(fullfile(chrpath,['chr' chrStr{c} '.fa']));
    ifa = imap(int16(fa.Sequence));        % numeric sequence
    Lchr(c) = numel(fa.Sequence);
    dn = find(diff(ifa)~=0);          % base transition positions
    lhp = diff(dn);              % hpol length between transitions
    ihp = ifa(dn(1:end-1)+1);
    for n = 1:4
        hpl(c,n,:) = hist(lhp(ihp==n),hpll);
    end
end
toc;
% Elapsed time is 181.168190 seconds.
%%
hpf = squeeze(nansum(hpl,1))/nansum(hpl(:));
hpf = hpf(:,2:end);
LL = hpll(2:end);
baseorder = 'ACGT';
save('HPfreq.mat','hpf','LL','baseorder','imap');

%%
figure(1); clf;
netHP = sum(squeeze(sum(hpl,1)));
footprintHP = netHP.*hpll;
cfootprintHP = cumsum(footprintHP,'reverse');
cfootprintHP(1) = NaN;
semilogy(hpll,cfootprintHP/sum(Lchr),'b-','linewidth',3);
% hold on;
% % semilogy(hpll,footprintHP/sum(Lchr),'r-');
% plot(hpll(11),cfootprintHP(11)/sum(Lchr),'bo','markersize',5,'linewidth',4)
% hold off;
set(gca,'xgrid','on','ygrid','on','fontsize',18,'xlim',[0 25]);
set(gca,'ytick',10.^(-4:0),'yticklabel',100*(1-10.^(-4:0)));
xlabel('Maximum Homopolymer Length');
ylabel('Fraction of Genome (%)');
title('Cumulative Fraction of hg38 Excluding Homopolymers');
%% HP as fraction of genome and flows
figure(2); clf;
netHP = squeeze(sum(hpl,1));
footprintHP = netHP.*repmat(hpll,4,1);
cfootprintHP = cumsum(footprintHP,2,'reverse')/sum(footprintHP(:));
cfootprintHP(:,1) = NaN;
subplot(1,2,1);
semilogy(hpll,cfootprintHP,'linewidth',1);
hold on;
semilogy(hpll,sum(cfootprintHP,1),'b-','linewidth',3);
hold off;
% hold on;
% % semilogy(hpll,footprintHP/sum(Lchr),'r-');
% plot(hpll(11),cfootprintHP(11)/sum(Lchr),'bo','markersize',5,'linewidth',4)
% hold off;
set(gca,'xgrid','on','ygrid','on','fontsize',18,'xlim',[0 25]);
set(gca,'ytick',10.^(-6:0),'yticklabel',100*10.^(-6:0),'ylim',10.^[-6 0]);
xlabel('Minimum Homopolymer Length');
ylabel('Cumulative Fraction of hg38 (%)');
title('Homopolymers in Human Genome');
subplot(1,2,2);
hpk = squeeze(sum(hpl));
chpk = cumsum(hpk,2,'reverse')/sum(hpk(:));
chpk(:,1) = NaN;
semilogy(hpll,chpk','linewidth',1);
hold on;
semilogy(hpll,sum(chpk),'b-','linewidth',3);
hold off;
set(gca,'xgrid','on','ygrid','on','fontsize',18,'xlim',[0 25]);
set(gca,'ytick',10.^(-6:0),'yticklabel',100*10.^(-6:0),'ylim',10.^[-6 0]);
legend({'A','C','G','T','All'},'location','NE');
xlabel('Minimum Homopolymer Length');
ylabel('Cumulative Fraction Flows (%)');
%%
figure(3); clf;
netHP = squeeze(sum(hpl,1));
footprintHP = netHP.*repmat(hpll,4,1);
cfootprintHP = cumsum(footprintHP,2,'reverse')/sum(footprintHP(:));
cfootprintHP(:,1) = NaN;
footprintHP = footprintHP/sum(footprintHP(:));
subplot(1,2,1);
semilogy(hpll,footprintHP,'linewidth',1);
hold on;
semilogy(hpll,sum(footprintHP,1),'b-','linewidth',3);
hold off;
% hold on;
% % semilogy(hpll,footprintHP/sum(Lchr),'r-');
% plot(hpll(11),cfootprintHP(11)/sum(Lchr),'bo','markersize',5,'linewidth',4)
% hold off;
set(gca,'xgrid','on','ygrid','on','fontsize',18,'xlim',[0 25]);
% set(gca,'ytick',10.^(-6:0),'yticklabel',100*10.^(-6:0),'ylim',10.^[-6 0]);
xlabel('Homopolymer Length');
ylabel('Fraction of hg38');
title('Homopolymers in Human Genome');
subplot(1,2,2);
hpk = squeeze(sum(hpl));
chpk = cumsum(hpk,2,'reverse')/sum(hpk(:));
chpk(:,1) = NaN;
hpk = hpk / sum(hpk(:));
semilogy(hpll,hpk','linewidth',1);
hold on;
semilogy(hpll,sum(hpk),'b-','linewidth',3);
hold off;
set(gca,'xgrid','on','ygrid','on','fontsize',18,'xlim',[0 25]);
% set(gca,'ytick',10.^(-6:0),'yticklabel',100*10.^(-6:0),'ylim',10.^[-6 0]);
title('Homopolymers in Human Genome');
legend({'A','C','G','T','All'},'location','NE');
xlabel('Homopolymer Length');
ylabel('Fraction Flows');

%%
hpc = squeeze(sum(hpl,1))';
fprintf('%6s %12s %12s %12s %12s\n', 'Length','A','C','G','T');
for i = 2:size(hpc,1)
    fprintf('%6d %12s %12s %12s %12s\n', hpll(i), ...
        num2cstr(hpc(i,1)), num2cstr(hpc(i,2)), num2cstr(hpc(i,3)), num2cstr(hpc(i,4)));
end

