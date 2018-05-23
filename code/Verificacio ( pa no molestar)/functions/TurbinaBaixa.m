function [ T,TAU,PI,P ] = TurbinaBaixa(T,alpha,CP,f,ETA,gam,P,PI,TAU,isTurboProp )
if isTurboProp
 TAU.tL = 1/(TAU.cH*TAU.r)*(TAU.r*TAU.cH*(TAU.r-1)+TAU.lamb)/(TAU.lamb-TAU.r*TAU.c+TAU.r);
 T.t5=TAU.tL*T.t45;
else
T.t5=T.t45-(T.t25-T.t2)*(((1+alpha)*CP.cold)/((1+f)*CP.hot*ETA.mec));
TAU.tL=T.t5/T.t45;
end
Ti.t5=T.t45-(T.t45-T.t5)/ETA.tL;
TAUi.tL=Ti.t5/T.t45;
PI.tL=tau2pi(TAUi.tL,gam.hot);
P.t5=P.t45*PI.tL;
PI.t=PI.tL*PI.tH;
TAU.t=TAU.tH*TAU.tL;
end

