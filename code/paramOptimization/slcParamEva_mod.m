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
ETA.f=0.98; %Suposici�
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

%% Considerem situaci� de Optimum Bypass ratio. Minimitza consum especific. (pag 370)
inc_c=0.4; %Increment del valor de PI.c
max_c = 30; %M�xim PI.c assumible
inc_f=0.1; 
max_f=6;
PI.c = 1 : inc_c : max_c; %Vector amb els valors de PI.c
PI.f= 1 : inc_f : max_f;
% i=1;
%C�lcul de alpha optima (pag 300 libro): 
% for pc = PI.c
%     TAU.c(i) = pi2tau(pc, gam.cold); 
%     for pf=PI.f
%        TAU.f(i) = pi2tau(pf, gam.cold);
%        alpha(i,j)=1/(TAU.r*TAU.f(i))*(TAU.lamb-TAU.r*(TAU.c(i)-1)-TAU.lamb/(TAU.r*TAU.c(i))-1/4*(sqrt(TAU.r*TAU.f(i)-1)+sqrt(TAU.r-1))^2);
%        Fadim(i,j)=a0/gc*(1+2*alpha(i,j))/(2*(1+alpha(i,j)))*(sqrt(2/(gam.cold-1)*(TAU.r*TAU.f(i)-1))-M0);
%         j = j + 1;
%     end
%     j = 1;
%     i = i + 1;
% end
%% PLOTS
% figure();
% surf(PI.c,PI.f,Fadim');
% xlabel('\pi_c'); ylabel('\pi_f'); zlabel('$$\hat{F}$$','Interpreter','Latex');
% figure();
% surf(PI.c,PI.f,alpha');
% xlabel('\pi_c'); ylabel('\pi_f'); zlabel('\alpha');
%% pc_opt and pf_opt iteration (search optimum values)
%initial values
pf_i(1) = PI.f(end);
pc_i(1) = 0;

j = 2; % for scrolling all pf_i values 
k = 2; % for scrolling all pc_i values 
diff_pf = 1; % difference between pf_opt in each iteration
TOL_pf = 0.2;% Maximum diff_pf value acceptable
diff_pc = 1; % difference between pf_opt in each iteration
TOL_pc = 0.2;% Maximum diff_pf value acceptable

%{
   Iteration stops when differnece between iterations of pc_opt and pf_opt
   is less than tolerance accepted.
%}
while  abs(diff_pf) > TOL_pf && abs(diff_pc) > TOL_pc 
    i=1;
    Fadim_anterior = 0; %always 0 for tarting
    %look for pc_opt
    for pc = PI.c
        TAU.c(i) = pi2tau(pc, gam.cold);
        pf = pf_i(j-1);
        TAU.f(i) = pi2tau(pf, gam.cold);
        alpha(i)=1/(TAU.r*TAU.f(i))*(TAU.lamb-TAU.r*(TAU.c(i)-1)-TAU.lamb/(TAU.r*TAU.c(i))-1/4*(sqrt(TAU.r*TAU.f(i)-1)+sqrt(TAU.r-1))^2);
        Fadim(i)=a0/gc*(1+2*alpha(i))/(2*(1+alpha(i)))*(sqrt(2/(gam.cold-1)*(TAU.r*TAU.f(i)-1))-M0);
        if (Fadim(i)-Fadim_anterior)/Fadim(i) < 0.0001
            pc_opt = pc; %optimum found
            %save convergence values for later comparation
            pc_i(k) = pc_opt;
            diff_pc = pc_i(k) - pc_i(k-1);
            k = k + 1;
            break
        else
            Fadim_anterior = Fadim(i);
            pc_opt = 0;
            pc_i(k) = pc_opt;
        end
        i = i + 1;
    end
    
    pf_opt = pf_i(j-1);
    fields = {'c','f'};
    clear alpha; clear Fadim;
    TAU = rmfield(TAU, fields);
    i = 1;
    for pf = PI.f
        pcf = pc_opt;
        TAU.c(i) = pi2tau(pcf, gam.cold);
        TAU.f(i) = pi2tau(pf, gam.cold);
        alpha(i)=1/(TAU.r*TAU.f(i))*(TAU.lamb-TAU.r*(TAU.c(i)-1)-TAU.lamb/(TAU.r*TAU.c(i))-1/4*(sqrt(TAU.r*TAU.f(i)-1)+sqrt(TAU.r-1))^2);
        Fadim(i)=a0/gc*(1+2*alpha(i))/(2*(1+alpha(i)))*(sqrt(2/(gam.cold-1)*(TAU.r*TAU.f(i)-1))-M0);
        if abs((Fadim(i)-Fadim_anterior)/Fadim(i)) < 0.05
            pf_opt = pf;
            break
        else
            Fadim_anterior = Fadim(i);
        end
        i = i + 1;
    end
    pf_i(j) = pf_opt;
    diff_pf = pf_i(j) - pf_i(j-1);
    j = j + 1;
end

fprintf('The optimum values are: \n pi_f = %.2f\n pi_c = %.2f\n alpha = %.2f\n', pf_opt, pc_opt,alpha(i));
%PI.f = pf_opt; PI.c = pc_opt; alpha = alpha(i);
