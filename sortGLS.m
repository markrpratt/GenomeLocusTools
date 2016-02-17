function [sR,is] = sortGLS_(inR)
sR = inR;
[sR.R,is] = sortrows(inR.R);
return;
