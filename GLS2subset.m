function sL = GLS2subset(L,ii)
sL = L;
sL.R = sL.R(ii,:);
return;
