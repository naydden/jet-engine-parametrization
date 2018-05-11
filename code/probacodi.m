%Per provar el codi de realengine.
PI.f=1;
PI.c=25;
T.t4=1700;
alpha=1;
T0=245;
P0=60000;
M0=0.8;
gam.cold=1.4;
gam.hot=1.4;
TAU.gam=6.93;
TAU.r=1.128;
PI.r=pi2tau(TAU.r,gam.cold);
T.t0=276.36;
P.t0=912000;
CP.hot=1005;
CP.cold=1005;
h=43e6; 
PI.d=0.96;
ETA.f=0.98; %Suposició
ETA.cH=0.98;
PI.b=0.94;
ETA.b=0.99;
ETA.tH=0.87;
ETA.tL=0.87;
PI.n=0.98;
ETA.mec=0.99;
[ PI,TAU,T,P,f ] =realengine( PI,TAU,T0,P0,M0,T,P,CP,gam,alpha,ETA,h);
mixer=false;
if mixer == true
  %COMPUTE THE MIXER  
else
    P.t6=P.t5;
    T.t6=T.t5;
    P.t16=P.t13;
    T.t16=P.t13;
end
[ P9,T9,P19,T19,Fadim,M9,M19 ] = Nozzle( T,P,PI,ETA,P0,T0,alpha,gam,f,M0);
