function [ T,TAU,PI,P ] = TurbinaBaixa(T,alpha,CP,f,ETA,gam,P,PI,TAU )
T.t5=T.t45-(T.t25-T.t2)*(((1+alpha)*CP.cold)/((1+f)*CP.hot*ETA.mec));
TAU.tL=T.t5/T.t45;
Ti.t5=T.t45-(T.t45-T.t5)/ETA.tL;
TAUi.tL=Ti.t5/T.t45;
PI.tL=tau2pi(TAUi.tL,gam.hot);
P.t5=P.t45*PI.tL;
PI.t=PI.tL*PI.tH;
TAU.t=TAU.tH*TAU.tL;
end

