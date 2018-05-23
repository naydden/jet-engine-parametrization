function [ T,TAU,PI,P ] = TurbinaBaixaTP(T,alpha,CP,f,ETA,gam,P,PI,TAU, T0)
theta0 = T.t0/T0;
thetat = T.t4/T0;
TAU.tL = 1/(TAU.c*theta0)*(theta0*TAU.c*(theta0-1)+thetat)/(thetat+theta0-theta0*TAU.c);
PI.tL=tau2pi(TAU.tL,gam.hot);
% P.t5=P.t45*PI.tL;
PI.t=PI.tL*PI.tH;
TAU.t=TAU.tH*TAU.tL;
T.t5 = TAU.tL*T.t45;

end

