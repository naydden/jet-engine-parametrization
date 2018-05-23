%% PARAMETRES GENERALS DONATS
%Thrusts
F=50000; %N

%Velocitat
v0 = 0; %[m/s] 600

%Camara combustio
T.t4=1143;
h=43e6; 

%Condicions ambientals
T0=288.15; %[K]
P0=101325; %[Pa]
rho0=1.225; %[kg/m^3]
Rgas=287;
% a0=sqrt(rho0*Rgas*T0); %301.6361
a0 = 340.26;%[m/s^2]
%     %atmosisa
%     height = 9500; %[m]
%     [T0, a0, P0, rho0] = atmosisa(height);

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
M0 = v0/a0;
T.t0 = T0*(1+(gam.cold-1)/2*M0^2);
P.t0=P0*(1+(gam.cold-1)/2*M0^2)^(gam.cold/(gam.cold-1));
TAU.lamb = T.t4/T0; 
%---------------------
TAU.r = T.t0/T0;
PI.r=tau2pi(TAU.r,gam.cold);

