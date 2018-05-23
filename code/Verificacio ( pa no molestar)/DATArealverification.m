%% PARAMETRES GENERALS DONATS
%Camara combustio
T.t4=1700;
h=43e6; 

%Condicions ambientals
T0=245; %[K]
P0=6e+04; %[Pa]

gam.cold=1.4;
CP.cold=1004;
R.cold = CP.cold-CP.cold/gam.cold;

gam.hot=1.4;
CP.hot=1004;
R.hot = CP.hot-CP.hot/gam.hot;

gc = 1;
%In engineering, gc is a unit conversion factor used to convert mass to force or vice versa

%Parametres del motor real
PI.d=0.92;
ETA.f=1; %Suposicio
ETA.cH=0.86;
PI.b=0.97;
ETA.b=0.98;
ETA.tH=0.92;
ETA.tL=1;
PI.n=0.92;
ETA.mec=1;

%Calculations
M0 = 0.8;
T.t0 = T0*(1+(gam.cold-1)/2*M0^2);
P.t0=P0*(1+(gam.cold-1)/2*M0^2)^(gam.cold/(gam.cold-1));
TAU.lamb = T.t4/T0; 
%---------------------
TAU.gam=T.t4/T0; %OJO, �s igual a TAU.lamb, nom�s canvia el nom
%---------------------
TAU.r = T.t0/T0;
PI.r=tau2pi(TAU.r,gam.cold);

