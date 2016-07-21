function density = densityGLS_(R, W, weight,mode)
%% density = densityGLS_(R, W, weight,mode)
% compute the density of GLS features in neighborhood of each R (start:stop) record
% R is the input segment start:stop matrix (nx2)
% W is the window width
% weight is the weight for each record
% mode controls behavior of counting.  mode=0 places all weight at start
% position of range, mode=1 places weight on start:stop

density = zeros(size(R,1),1,'single');
mnR = min(R(:,1))-W;
mxR = max(R(:,2))+W;
lR = mxR-mnR+1;
D = zeros(lR,1,'single');

% apply weights to segment vector
switch mode
    case 0
        for i = 1:size(R,1)
            cc = R(i,1) - mnR;
            D(cc) = D(cc) + weight(i);
        end;
    case 1
        for i = 1:size(R,1)
            cc = (R(i,1):R(i,2)) - mnR;
            D(cc) = D(cc) + weight(i);
        end;
end;

% average weights in window for density
ww = int32(-round(W/2):round(W/2));
for i = 1:size(R,1)
    cc = int32(round(mean(R(i,1:2)))) + ww - mnR;  % window about center of range for now
    density(i) = mean(D(cc));
end;
        
return;