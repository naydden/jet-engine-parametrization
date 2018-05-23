function [ P,f ] = CambraCombustioTP( PI,P,CP,TAU,gam,ETA,h,T0, T)
%Apartat c exercici 33

P.t4=PI.b*P.t3; %per definició

%parametres adim
theta0 = T.t0/T0;
thetat = T.t4/T0;

f = CP.hot*T0/h*(thetat-theta0*TAU.c);

end

