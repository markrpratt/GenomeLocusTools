function cstrs = GLS2igvstr(L,ii,verbose)
if nargin < 2
    ii = 1:size(L.R,1);
end;
ii = ii(:)';
if nargin < 3
    verbose = true;
end;
cstr = cell(size(ii));
for i = ii
    cstr{i} = sprintf('chr%s:%d-%d', L.segNames{L.R(i,1)}, L.R(i,2:3));
    if verbose
        fprintf('%s\n', cstr{i});
    end;
end;
