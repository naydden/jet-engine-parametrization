function [ m0,mf ] = Fluxosmasics( f,Fadim,F,a0)
m0=F/(Fadim*a0);
mf=f*m0;


end

