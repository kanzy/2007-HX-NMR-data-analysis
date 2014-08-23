%%hxfit.m 2007-07-10: called by normvol.m to fit HX rate from gCOSY NMR data (vol measured in Felix):  
%%current version updated at 2007-07-21 

function f = hxfit(x, time, i, intenCorr)

A=x(1);
K=x(2);

funxdata = A*exp(-K*time*3600);

ydata = intenCorr(i,:);
     
   f = ydata - funxdata;