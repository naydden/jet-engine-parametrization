function [ m0,mf,msec ] = Fluxosmasics( f,Fadim,F,a0,alpha)
m0=F/(Fadim*a0); %Flux masic del core
mf=f*m0; %Flux masic del combustible
msec=m0*alpha; %Flux secundari

end

