function [ PI,P,TAU,T ] = Compressor( PI,P,gam,T,TAU,ETA )
PI.cH=PI.c/PI.f;
TAU.c=pi2tau(PI.c,gam.cold);
P.t3=P.t25*PI.cH;
TAUi.cH=pi2tau(PI.cH,gam.cold);
Ti.t3=T.t25*TAUi.cH;
T.t3=T.t25+(Ti.t3-T.t25)/ETA.cH;
TAU.cH=T.t3/T.t25;
end

