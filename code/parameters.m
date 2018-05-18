clear; clc; close all;
%{
%}
%C = @tau2pi;
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


%% OPTIMIZATION
%input {PI.c, alpha}
inc_c = 0.4;
inc_a = 0.1;

max_alpha = 5;
max_c = 30;

i = 1;
j = 1;
PI.c = 1 : inc_c : max_c;

for pc = PI.c
    TAU.c(i) = pi2tau(pc, gam.cold); 
    for alpha = 0: inc_a : max_alpha
        TAU.f(i,j) = (TAU.lamb - TAU.r*(TAU.c(i) - 1) - TAU.lamb/(TAU.r*TAU.c(i)) + alpha*TAU.r+1)*1/(TAU.r*(1+alpha));
        PI.f(i,j) = tau2pi(TAU.f(i,j), gam.cold);
        j = j + 1;
    end
    j = 1;
    i = i + 1;
end
%alpha = 0: inc_a : max_alpha;
[I,J] = size(PI.f);
%% PLOT
% surf(PI.c, alpha, PI.f')
% xlabel('\pi_c'); ylabel('\alpha'); zlabel('\pi_f');

%% Adimensional F vs. PI.c
%book page 461
R.c = (gam.cold-1)/gam.cold*CP.cold;
R.t = (gam.hot-1)/gam.hot*CP.hot;
V.v0 = a0*M0;
if M0 <= 1
    ETA.r = 1;
else
    ETA.r = 1-0.075*(M0-1)^1.35;
end
i = 1;
j = 1;
f = zeros(size(PI.c));
M9 = zeros(I,J);
M19 = zeros(I,J);
gc = 6.674e-11;%[ m3?kg?1?s?2] constante newton

for pc = PI.c
    f(i) = (TAU.lamb - TAU.r*TAU.c(i))/(ETA.b*h/(CP.cold*T0)-TAU.lamb);
    for alpha = 0: inc_a : max_alpha
        TAU.t(i,j) = 1-1/(ETA.mec*(1+f(i)))*TAU.r/TAU.lamb*(TAU.c(i)-1+alpha*(TAU.f(i,j)-1));
        PI.t(i,j) = tau2pi(TAU.t(i,j), gam.hot);
        
        %core
        P.P0_P9(i,j) = 1;
        %compte, dona <1 pero no té sentit
        P.t9_9(i,j) = P.P0_P9(i,j)*PI.r*PI.d*PI.c(i)*PI.b*PI.t(i,j)*PI.n ; %Suposem adaptada P0/P9 = 1
%         if P.t9_9(i,j) < 1
%             warning('%d, %d',i,j)
%         end
        M9(i,j) = sqrt(2/(gam.hot-1)*(P.t9_9(i,j)^((gam.hot-1)/gam.hot)-1));
        if real(M9(i,j)) > 1 %tobera xocada
            M9(i,j) = 1;
            P.t9_9(i,j) = (1 + M9(i,j)^2*(gam.hot-1)/2)^(gam.hot/(gam.hot-1));    
        end
        % [rows,cols,vals] =find(M9(imag(M9)== 0))
        T.T9_T0(i,j) = TAU.lamb*TAU.t(i,j)*CP.cold/(pi2tau(P.t9_9(i,j),gam.hot)*CP.hot);
        V.v9_a0(i,j) = M9(i,j)*sqrt((gam.hot*R.t)/(gam.cold*R.c)*T.T9_T0(i,j));
        
        %secondary
        P.P0_P19(i,j) = 1;
        P.t19_19(i,j) = P.P0_P19(i,j)*PI.r*PI.d*PI.f(i,j)*PI.n ; %Suposem adaptada P0/P19 = 1
        M19(i,j) = sqrt(2/(gam.cold-1)*(P.t19_19(i,j)^((gam.cold-1)/gam.cold)-1));
        if real(M19(i,j)) > 1 %tobera xocada
            M19(i,j) = 1;
            P.t19_19(i,j) = (1 + M19(i,j)^2*(gam.cold-1)/2)^(gam.cold/(gam.cold-1));    
        end
        T.T19_T0(i,j) = TAU.lamb*TAU.f(i,j)/(pi2tau(P.t19_19(i,j),gam.cold));
        V.v19_a0(i,j) = M19(i,j)*sqrt(T.T19_T0(i,j));
        
        F.adim(i,j) = a0/((1+alpha)*gc)*((1+f(i))*V.v9_a0(i,j)-M0+(1+f(i))*...
            (R.t*T.T9_T0(i,j)*(1-P.P0_P9(i,j)))/(R.c*V.v9_a0(i,j)*gam.cold))+...
            alpha*a0/((1+alpha)*gc)*(V.v19_a0(i,j)-M0+...
            T.T19_T0(i,j)/V.v19_a0(i,j)*(1-P.P0_P19(i,j))/gam.cold);  
        j = j + 1;
    end
    j = 1;
    i = i + 1;
end
%% Process F.adim
%look for real values
for i = 1:I
    for j = 1:J
        if isreal(F.adim(i,j))
            F.adim_c(i,j) = F.adim(i,j);
        else
            F.adim_c(i,j) = 0;
        end
    end
end
%{
j = 1;
TOL.val = 2.5/100;
alpha = 0: inc_a : max_alpha;

for i = 1:I
    TOL.c = 100;
    while TOL.c > TOL.val 
        if j <= J - 1
            TOL.c = ((F.adim_c(i,j+1)-F.adim_c(i,j))/...
                (alpha(j+1)-alpha(j)));
        else
            TOL.c = 0;
        end
        j = j + 1;
    end
    F.alpha(i) = alpha(j);
    j=1;
end
[row, column, value] = find(alpha == F.alpha(10));
% PI.c
maxPIc = find(F.adim_c(:,column) == max(F.adim_c(:,column)));
PI.c_slct = PI.c(maxPIc);
PI.f_slct = PI.f(maxPIc, column);
%}

%% Select max's line
for i = 1:I
    F.maxs(i) = max(F.adim_c(i,:));
 
    vol = find(F.adim_c(i,:) == max(F.adim_c(i,:)));
    F.alpha_pos(i) = vol(1);
end
i = 1;
TOL.val = 2.5/100;
alpha = 0: inc_a : max_alpha;
F.alpha = 0: inc_a : max_alpha;

TOL.c = 100;
while TOL.c > TOL.val
    if i <= I - 1
        TOL.c = ((F.maxs(i+1)-F.maxs(i))/...
            (PI.c(i+1)-PI.c(i)));
    else
        TOL.c = 0;
    end 
    F.PIc = (PI.c(i+1)+PI.c(i))/2;
    F.a = F.alpha(F.alpha_pos(i));
    i = i + 1;
end

F.PIf = PI.f(i-1, F.alpha_pos(i));

%% Plot
% figure()
% plot(alpha, F.adim_c');
% figure()
% plot(PI.c, F.adim_c);
figure();
surf(PI.c,alpha,F.adim_c'/a0)
xlabel('\pi_c'); ylabel('\alpha'); zlabel('\bar F');
figure();
plot(PI.c, F.maxs)
xlabel('\pi_c'); ylabel('$$\hat{F}$$','Interpreter','Latex');


