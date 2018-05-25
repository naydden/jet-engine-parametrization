clear; clc; close all;
%{
%}
%% Values
T.t4 = 1780; %[K]
height = 9500; %[m]
[T0, a0, P0, rho0] = atmosisa(height);
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


%% Design parameters
%vector creation
inc_c = 0.4;
inc_a = 0.1;
inc_f = 0.5;

max_alpha = 10;
max_c = 30;
max_f = 10;

PI.c = 1 : inc_c : max_c;
alpha = 0: inc_a : max_alpha;
PI.f = 0 : inc_f : max_f;

I = max(size(PI.c));
J = max(size(alpha));
K = max(size(PI.f));

%% Adimensional F vs. PI.c
%book page 347
R = (gam.cold-1)/gam.cold*CP.cold;
i = 1;
for pc = PI.c
    TAU.c(i) = pi2tau(pc, gam.cold);
    i = i + 1;
end
k = 1;
for pf = PI.f
    TAU.f(k) = pi2tau(pf, gam.cold);
    k = k + 1;
end

i = 1; j = 1; k = 1;
f = zeros(I,J,K);
S = zeros(I,J,K);

gc = 1;%gc engineering

for pc = PI.c
    for a = alpha
        for pf = PI.f
            v.v9_a0(i,j,k) = sqrt(2/(gam.cold -1)*(TAU.lamb-TAU.r*(pi2tau(pc,gam.cold))...
                - 1 + a*(pi2tau(pf,gam.cold)))-TAU.lamb/(TAU.r*pi2tau(pc,gam.cold)));
            v.v19_a0(i,j,k) = sqrt(2/(gam.cold -1)*(TAU.r*pi2tau(pf,gam.cold)-1));
            F.adim(i,j,k) = a0/gc*1/(1 + a)*(v.v9_a0(i,j,k) - M0 + a*...
                (v.v19_a0(i,j,k) - M0));
            f(i,j,k) = CP.cold*T0/h*(TAU.lamb-TAU.r*pi2tau(pc,gam.cold));
            S(i,j,k) = f(i,j,k)/((1+a)*(F.adim(i,j,k)));
            k = k + 1;
        end
        k = 1;
        j = j + 1;
    end
    k = 1;
    j = 1;
    i = i + 1;
end
%% Complex erasing
%look for real values
for i = 1:I
    for j = 1:J
        for k = 1:K
            if isreal(F.adim(i,j,k))
                F.adim_c(i,j,k) = F.adim(i,j,k);
            else
                F.adim_c(i,j,k) = 0;
            end
        end
    end
end

%% Plots
for k = 1:K
    hold on
    plot(PI.f, F.adim_c(:,40,k));
end
hold off


% plot(PI.c, F.adim_c(:,40,10)/a0);


% for k = 1:K
%     kPlot(k) = F.adim_c(50,40,k);
% end
% plot(PI.f, kPlot);