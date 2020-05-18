%% Tabulate homopolymers in human genome
clear;
chrpath = 'C:\Users\MarkPratt\Documents\Genome\hg38\chroms';
chrStr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13', ...
    '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M'};
imap(int16('ACGTacgtNn')) = int8([0 1 2 3 0 1 2 3 7 7]);
cmap = 'ACGT';              % need to add one to above index
k3 = int16([1 4 16]);
hh = 0:63;
H = zeros(numel(hh),numel(chrStr));   % histogram of 3-mer frequency by Chr
Lchr = zeros(size(chrStr));
tic;
parfor c = 1:numel(chrStr)
    fa = fastaread(fullfile(chrpath,['chr' chrStr{c} '.fa']));
    Lchr(c) = numel(fa.Sequence);
    ifa = imap(int16(fa.Sequence));        % numeric sequence
    C = conv2(ifa,k3,'same');
    iok = ifa<7;
    iok = imerode(iok,[1 1 1]);
    H(:,c) = hist(C(iok),hh);
end
toc;
% Elapsed time is 216.778428 seconds.
%% All possible 3-mers
nH = 64*H ./ repmat(sum(H),numel(hh),1);
Tmer = cell(64,1);
for i = 1:64
    j = int16(i-1);
    Tmer{i} = cmap([mod(j,4) mod(bitshift(j,-2),4) bitshift(j,-4)]+1);
end
gH = 64*sum(H,2)/sum(H(:));
[gH,is] = sort(gH,'ascend');       % sort by genome average
sH = nH(is,:);
sTmer = Tmer(is);
disp(Tmer);

figure(1); clf;
cm = jet(numel(chrStr));
for i = 1:numel(chrStr)
    plot(sH(:,i),'-.','color',cm(i,:),'linewidth',1);
    hold on;
end
plot(gH,'-','linewidth',2,'color',[1 1 1]*0.35);
hold off;
colormap hsv;
set(gca,'xtick',hh+1,'xticklabel',sTmer,'xgrid','on','ygrid','on','xlim',[1 64],'fontsize',12);
xtickangle(45);
legStr = chrStr;
legStr{end+1} = 'Genome';
legend(legStr,'location','nw');
title('Distribution of 3-mers in hg38');
ylabel('Relative Frequency');

for i = 1:64
    fprintf('%s %0.4f\n', sTmer{i}, gH(i));
end