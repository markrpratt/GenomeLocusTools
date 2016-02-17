function L = binary2GLS_(V)
%% L = binary2GLS_(V)
% returns Nx2 array of start/stop positions for true values in V

d = diff(V);
ii = find(d);
ii(d(ii)>0) = ii(d(ii)>0)+1;

if V(1) == 1
    if V(end) == 1
        ii = [1 ii numel(V)];
    else
        ii = [1 ii];
    end;
else
    if V(end) == 1
        ii = [ii numel(V)];
    end;
end;

L = int32(reshape(ii,2,numel(ii)/2)');
   