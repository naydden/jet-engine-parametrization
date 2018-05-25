function [ T,TAU,PI,P ] = TurbinaAltaTP( T,CP,ETA,f,gam,TAU,P,PI)
T.t45=T.t4-(T.t3-T.t2)*(CP.cold/(ETA.mec*(1+f)*CP.hot));
TAU.tH=T.t45/T.t4;
Ti.t45=T.t4-(T.t4-T.t45)/ETA.tH;
TAUi.tH=Ti.t45/T.t4;
PI.tH=tau2pi(TAUi.tH,gam.hot);
P.t45=P.t4*PI.tH;
end

