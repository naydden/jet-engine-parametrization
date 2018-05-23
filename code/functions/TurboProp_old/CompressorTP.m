function [ PI,P,TAU,T ] = CompressorTP( PI,P,gam,T,TAU,ETA )
%IN: PI.c, gam, T.t2, ETA.cH, ETA.f
%OUT:TAU.c, T.t3, P.t3

TAU.c=pi2tau(PI.c,gam.cold);
%Per definició - No s'hi ha posat rendiment
P.t3=P.t2*PI.c;
T.t3=T.t2*TAU.c;
end
