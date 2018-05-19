function [ P,f ] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0)
P.t4=PI.b*P.t3;
f=(CP.hot*TAU.gam-CP.cold*TAU.r*TAU.d*TAU.cH*TAU.f)/((ETA.b*h/T0)-CP.hot*TAU.gam);
end

