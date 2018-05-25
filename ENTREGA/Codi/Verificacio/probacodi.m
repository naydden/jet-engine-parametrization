%Per provar el codi de realengine.
F=25000; %N
PI.f=1.4;
PI.c=25;
T.t4=1700;
alpha=1;
T0=245;
P0=60000;
M0=0.8;
rho0=1.225;
Rgas=287;
a0=sqrt(rho0*Rgas*T0);
gam.cold=1.4;
gam.hot=1.4;
TAU.r=1.128;
TAU.gam=T.t4/T0;
PI.r=pi2tau(TAU.r,gam.cold);
T.t0=276.36;
P.t0=912000;
CP.hot=1005;
CP.cold=1005;
h=43e6; 
%Parámetres del motor real
PI.d=0.96;
ETA.f=0.98; %Suposició
ETA.cH=0.98;
PI.b=0.94;
ETA.b=0.99;
ETA.tH=0.87;
ETA.tL=0.87;
PI.n=0.98;
ETA.mec=0.99;
%Cálcul de les etapes del jet:
[T,P,TAU] = Difusor( T,P,TAU,PI );
[P,TAU,T] = Fan( P,PI,gam,ETA,T,TAU );
[PI,P,TAU,T] = Compressor( PI,P,gam,T,TAU,ETA);
[P,f] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0 );
[T,TAU,PI,P] = TurbinaAlta( T,CP,ETA,f,gam,TAU,P,PI);
[ T,TAU,PI,P ] = TurbinaBaixa(T,alpha,CP,f,ETA,gam,P,PI,TAU);
mixer=false;
%Punts d'entrada a mixer: 1.3 i 5. Punt a la sortida del mixer: 6
if mixer == true
  %COMPUTE THE MIXER  
  [ T,P,M9 ] = ToveraPrimari( P0,T,PI,gam,P,TAU,mixer);
else
    P.t6=P.t5;
    T.t6=T.t5;
    P.t16=P.t13;
    T.t16=P.t13;
    TAU.b=tau2pi(PI.b,gam.hot);
    TAU.n=tau2pi(PI.n,gam.hot);
    [ T,P,M9 ] = ToveraPrimari( P0,T,PI,gam,P,TAU,mixer );
    [ T,P,M19] = Toverasecundari( T,PI,P,P0,gam,TAU );
end
Fadim = Fadimensional( f,M9,M19,alpha,T,P,gam,mixer,T0,M0,P0);
[ m0,mf ] = Fluxosmasics( f,Fadim,F,a0);
%Cálcul árees:
M5=1;
m5=m0+mf;
A.e5 = Area(M5,gam.hot,P.t5,T.t5,Rgas,m5);
