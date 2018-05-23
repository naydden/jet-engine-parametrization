function [ PI,P,TAU,T ] = Compressor( PI,P,gam,T,TAU,ETA,isTurboProp )
%IN: PI.c, PI.f, PI.cH, gam, T.t25, ETA.cH
%OUT:TAU.cH, T.t3, P.t3, TAU.c, PI.cH
PI.cH=PI.c/PI.f;
if isTurboProp %No hi ha fan així que els valors 25 son iguals als valors 2
    P.t25=P.t2;
    T.t25=T.t2;
    PI.cH=PI.c;
end
TAU.c=pi2tau(PI.c,gam.cold);
P.t3=P.t25*PI.cH;
TAUi.cH=pi2tau(PI.cH,gam.cold);
Ti.t3=T.t25*TAUi.cH;
T.t3=T.t25+(Ti.t3-T.t25)/ETA.cH;
TAU.cH=T.t3/T.t25;
end

