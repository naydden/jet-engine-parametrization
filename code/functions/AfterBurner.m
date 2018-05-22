function [ P,f_ab ] = AfterBurner( PI,P,CP,TAU,gam,ETA,h,T0,R)

gam.AB = gam.hot;
CP.AB = CP.hot;

R.AB = (gam.AB - 1)/gam.AB*CP.AB;



P.t4=PI.b*P.t3;
f_ab=(CP.hot*TAU.gam-CP.cold*TAU.r*TAU.d*TAU.cH*TAU.f)/((ETA.b*h/T0)-CP.hot*TAU.gam);

end

