clear; clc; close all;
%{
%}
%C = @tau2pi;
%% Values
T.t4 = 1780; %[K]
h = 9500; %[m]
[T0, a0, P0, rho0] = atmosisa(h);
v0 = 600/3.6; %[m/s]
gam.cold = 1.4;
gam.hot=1.3;
CP.hot=1200;
CP.cold=1004;

% Eficiencias i ratis
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

%Calculations
M0 = v0/a0;
T.t0 = T0*(1+(gam.cold-1)/2*M0^2);
P.t0=P0*(1+(gam.cold-1)/2*M0^2)^(gam.cold/(gam.cold-1));
TAU.lamb = T.t4/T0;
TAU.r = T.t0/T0;
PI.r=pi2tau(TAU.r,gam.cold);


%% OPTIMIZATION
%input {PI.c, alpha}
inc_c = 0.2;
inc_a = 0.2;

i = 1;
j = 1;
PI.c = 1 : inc_c : 60;

for pc = PI.c
    TAU.c(i) = pi2tau(pc, gam.cold); 
    for alpha = 0: inc_a : 30
        TAU.f(i,j) = (TAU.lamb - TAU.r*(TAU.c(i) - 1) - TAU.lamb/(TAU.r*TAU.c(i)) + alpha*TAU.r+1)*1/(TAU.r*(1+alpha));
        PI.f(i,j) = tau2pi(TAU.f(i,j), gam.cold);
        j = j + 1;
    end
    j = 1;
    i = i + 1;
end
alpha = 0: inc_a : 30;

%% PLOT
surf(PI.c, alpha, PI.f')
xlabel('\pi_c'); ylabel('\alpha'); zlabel('\pi_f');
