function L = Table2GLS(T)
%% Create genome locus structure and annotate 
L.refID = 'GRCh37';
L.hdr = [];
L.R = zeros(size(T,1),3,'int32');
[L.segNames, ~, L.R(:,1)] = unique(T.chrom);
% strip chr prefix, and change MT to M
L.segNames = regexprep(L.segNames,'chr','');
L.segNames = regexprep(L.segNames,'MT','M');
% standardize order of primary segments
StdSegNames = { ...
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', ...
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M'};
[~,ii,jj] = intersect(StdSegNames,L.segNames);
kk = setdiff(1:numel(L.segNames),jj);  % move non-std names in GLS to end
I = zeros(size(L.R,1),1,'int32');
for n = numel(kk):-1:1
   NewSegNames{n+numel(StdSegNames)} = L.segNames{kk(n)};
    I(L.L(:,1)==kk(n)) = n+numel(StdSegNames);
end;
% populate first 25 indices with standard chrs
NewSegNames(1:numel(StdSegNames)) = StdSegNames;
for n = 1:numel(jj)  
    I(L.R(:,1)==jj(n),1) = ii(n);
end;
L.segNames = NewSegNames;
L.R(:,1) = I;
L.R(:,2) = T.pos;
L.R(:,3) = T.pos;
