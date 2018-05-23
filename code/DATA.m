%% PAR�METRES GENERALS DONATS
%Thrusts
F=25000; %N

%Velocitat
v0 = 600/3.6; %[m/s] 600

%C�mara combusti�
T.t4=1780;
h=43e6; 

%Condicions ambientals
T0=226.4000; %[K]
P0=2.8524e+04; %[Pa]
rho0=0.4389; %[kg/m^3]
Rgas=287;
% a0=sqrt(rho0*Rgas*T0); %301.6361
a0 = 301.6361;%[m/s^2]
%     %atmosisa
%     height = 9500; %[m]
%     [T0, a0, P0, rho0] = atmosisa(height);

gam.cold=1.4;
CP.cold=1004;
R.cold = CP.cold-CP.cold/gam.cold;

gam.hot=1.3;
CP.hot=1200;
R.hot = CP.hot-CP.hot/gam.hot;

gc = 1;
%In engineering, gc is a unit conversion factor used to convert mass to force or vice versa

%Parametres del motor real
PI.d=0.96;
ETA.f=0.88; %Suposicio
ETA.cH=0.88;
PI.b=0.94;
ETA.b=0.99;
ETA.tH=0.87;
ETA.tL=0.87;
PI.n=0.98;
ETA.mec=0.99;

%Calculations
M0 = v0/a0;
T.t0 = T0*(1+(gam.cold-1)/2*M0^2);
P.t0=P0*(1+(gam.cold-1)/2*M0^2)^(gam.cold/(gam.cold-1));
TAU.lamb = T.t4/T0; 
%---------------------
TAU.gam=T.t4/T0; %OJO, �s igual a TAU.lamb, nom�s canvia el nom
%---------------------
TAU.r = T.t0/T0;
PI.r=tau2pi(TAU.r,gam.cold);
%% POST-COMBUSTOR
%After Burner parameters
gam.AB = gam.hot;
CP.AB = CP.hot;
T.t7 = 2400; %K
ETA.AB = 0.98;

%% TurboProp
ETA.prop = 0.7;
