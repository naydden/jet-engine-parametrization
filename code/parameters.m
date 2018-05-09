clear; clc; close all;
%{
%}
%C = @tau2pi;
%% Values
T.t4 = 1780; %[K]
h = 9500; %[m]
[T0, a0, P0, rho0] = atmosisa(h);
v0 = 600/3.6; %[m/s]
lamb.cold = 1.4;
lamb.hot=1.3;
CP.hot=1200;
CP.cold=1004;

%Calculations
M0 = v0/a0;
T.t0 = T0*(1+(lamb.cold-1)/2*M0^2);
P.t0=P0*(1+(lamb.cold-1)/2*M0^2)^(lamb.cold/(lamb.cold-1));
TAU.lamb = T.t4/T0;
TAU.r = T.t0/T0;


%% OPTIMIZATION
%input {PI.c, alpha}
inc_c = 0.2;
inc_a = 0.2;

i = 1;
j = 1;
PI.c = 1 : inc_c : 60;

for pc = PI.c
    TAU.c(i) = pi2tau(pc, lamb.cold); 
    for alpha = 0: inc_a : 30
        TAU.f(i,j) = (TAU.lamb - TAU.r*(TAU.c(i) - 1) - TAU.lamb/(TAU.r*TAU.c(i)) + alpha*TAU.r+1)*1/(TAU.r*(1+alpha));
        PI.f(i,j) = tau2pi(TAU.f(i,j), lamb.cold);
        j = j + 1;
    end
    j = 1;
    i = i + 1;
end
alpha = 0: inc_a : 30;

%% PLOT
surf(PI.c, alpha, PI.f')
xlabel('\pi_c'); ylabel('\alpha'); zlabel('\pi_f');
