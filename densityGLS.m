function density = densityGLS(GLS, W, weight,mode )
%% density = densityGLS(GLS, W, weight )
% compute the density of GLS features in neighborhood of each GLS record
% GLS is the input genome locus set
% W is the optional width, default is 100bp
% width is the optional weight for each record, default is 1
% mode controls behavior of counting.  mode=0 places all weight at start
% position of range, mode=1 places weight on start:stop
if nargin < 4
    mode = 0;
end;
if nargin < 3
    weight = ones(size(GLS.R,1),1,'single');
end;
if nargin < 2
    W = 100;
end;
density = zeros(size(GLS.R,1),1);
for i = 1:numel(GLS.segNames)
    ii = find(GLS.R(:,1)==i);
    density(ii) = densityGLS_(GLS.R(ii,2:3), W, weight(ii), mode);
end;
return;