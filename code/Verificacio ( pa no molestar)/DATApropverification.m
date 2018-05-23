%% PARAMETRES GENERALS DONATS

%Camara combustio

h=43e6; 

%Condicions ambientals
T0=249.15; %[K]
P0=4.7320e+04; %[Pa]
Rgas=287;

gam.cold=1.4;
CP.cold=1004;
R.cold = CP.cold-CP.cold/gam.cold;

gam.hot=1.4;
CP.hot=1004;
R.hot = CP.hot-CP.hot/gam.hot;

gc = 1;
%In engineering, gc is a unit conversion factor used to convert mass to force or vice versa

%Parametres del motor real
PI.d=1;
ETA.f=1; %Suposicio
ETA.cH=1;
PI.b=1;
ETA.b=1;
ETA.tH=1;
ETA.tL=1;
PI.n=1;
ETA.mec=1;

%Calculations
M0 = 0.5;
T.t0 = T0*(1+(gam.cold-1)/2*M0^2);
P.t0=P0*(1+(gam.cold-1)/2*M0^2)^(gam.cold/(gam.cold-1));
TAU.lamb = 7; 
T.t4=TAU.lamb*T0;
TAU.r = T.t0/T0;
PI.r=tau2pi(TAU.r,gam.cold);

%% TurboProp
ETA.prop = 95; %Rendiment turboprop
