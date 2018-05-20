clear; clc; close all;
%{
%}
%% Values
T.t4 = 1780; %[K]
height = 9500; %[m]
[T0, a0, P0, rho0] = atmosisa(height);
v0 = 600/3.6; %[m/s] 600
gam.cold = 1.4;
gam.hot=1.3;
CP.hot=1200;
CP.cold=1004;
gc = 1;
%In engineering, gc is a unit conversion factor used to convert mass to force or vice versa

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
inc_c = 0.4;
inc_a = 0.1;

min_alpha = 0;
max_alpha = 10;
max_c = 30;

i = 1;
j = 1;
PI.c = 1 : inc_c : max_c;

for pc = PI.c
    TAU.c(i) = pi2tau(pc, gam.cold); 
    for alpha = 0: inc_a : max_alpha
        TAU.f(i,j) = (TAU.lamb - TAU.r*(TAU.c(i) - 1) - TAU.lamb/(TAU.r*TAU.c(i)) + alpha*TAU.r+1)/(TAU.r*(1+alpha));
        PI.f(i,j) = tau2pi(TAU.f(i,j), gam.cold);
        V.v19_a0(i,j) = sqrt(2/(gam.cold-1)*(TAU.r*TAU.f(i,j)-1));
        F.adim(i,j) = a0/gc*(V.v19_a0(i,j)-M0);
        j = j + 1;
    end
    j = 1;
    i = i + 1;
end
alpha = 0: inc_a : max_alpha;
[I,J] = size(PI.f);

%% Optimum alpha
j = 1;
TOL.val = 2.5/100;
for i = 1:I
    TOL.c = 100;
    while TOL.c > TOL.val 
        if j <= J - 1
            TOL.c = ((F.adim(i,j)-F.adim(i,j+1)))/(F.adim(i,j));
        else
            TOL.c = 0;
            TOL.stat(i,j) = 0;
        end
        j = j + 1;
    end
    TOL.stat(i,j-1) = 1;
    F.alpha(i) = alpha(j-1);
    j=1;
end
%% Optimum PIc
i = 1;
TOL.val = 100/100;
for j = 1:J
    TOL.c = 100;
    while TOL.c > TOL.val 
        if i <= I - 1
            TOL.c = abs(((F.adim(i,j)-F.adim(i+1,j))));
        else
            TOL.c = 0;
            TOL.stat(i,j) = 0;
        end
        i = i + 1;
    end
    TOL.stat(i-1,j) = 1;
    F.PIc(i) = PI.c(i-1);
    i=1;
end
F.PIc = max(F.PIc);
F.a = F.alpha((PI.c == max(F.PIc)));
F.PIf = PI.f((PI.c == max(F.PIc)),(alpha == F.a));


%% PLOTs
figure();
surf(PI.c,alpha,F.adim');
xlabel('\pi_c'); ylabel('\alpha'); zlabel('$$\hat{F}$$','Interpreter','Latex');
figure();
surf(PI.c,alpha,PI.f');
xlabel('\pi_c'); ylabel('\alpha'); zlabel('\pi_f');
