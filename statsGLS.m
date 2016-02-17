function s = statsGLS(S,verbose)
%% function s = statsGLS(S, verbose)
% returns a 2 x n array where n is the number of segments in S
% rows are: number of ranges, number of bases, 
if nargin < 2
    verbose = false;
end;

if isempty(S)
    s = [0 0];
    return;
end;
s = zeros(numel(S.segNames),2);
for i = 1:numel(S.segNames)
    L = diff(S.R(S.R(:,1)==i,2:3),1,2)+1;
    s(i,1) = numel(L);
    s(i,2) = sum(L);
end;

if verbose
    strlen = max(7, max(cellfun(@length,S.segNames))+1);
    nr = sum(s(:,1));
    nb = sum(s(:,2));
    nrlen = max(7,1+floor(log10(nr)));
    nblen = max(6,1+floor(log10(nb)));
    fprintf('%*s  %*s  %*s\n', strlen, 'Segment', nrlen, 'nRanges', nblen, 'nBases');
    for i = 1:size(s,1)
        fprintf('%*s  %*d  %*d\n', strlen, S.segNames{i}, ...
            nrlen, s(i,1), nblen, s(i,2));
    end;
    fprintf('%*s  %*d  %*d\n', strlen, 'All', nrlen, nr, nblen, nb);
end;
return;
