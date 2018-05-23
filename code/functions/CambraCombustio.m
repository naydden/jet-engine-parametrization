function [ P,f ] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0,isTurboProp)
if isTurboProp %No hi ha fan
    TAU.f=1;
end
P.t4=PI.b*P.t3;
f=(CP.hot*TAU.lamb-CP.cold*TAU.r*TAU.d*TAU.cH*TAU.f)/((ETA.b*h/T0)-CP.hot*TAU.lamb);
end

