function [ PI, alpha ] = opt_parameters( M0, a0, gam, gc, PI, TAU )
%TURBOFAN optimal parameters
%   Compute PI.f, PI.c, alpha
%% Considerem situació de Optimum Bypass ratio. Minimitza consum especific. (pag 370)
inc_c=0.4; %Increment del valor de PI.c
max_c = 30; %Màxim PI.c assumible
inc_f=0.1;
max_f=6;
PI.c = 1 : inc_c : max_c; %Vector amb els valors de PI.c
PI.f= 1 : inc_f : max_f;

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
while  abs(diff_pf) > TOL_pf || abs(diff_pc) > TOL_pc
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
PI.f = pf_opt; PI.c = pc_opt; alpha = alpha(i);



end

