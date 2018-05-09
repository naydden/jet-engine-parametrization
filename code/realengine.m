function [ PI,TAU,T,P,f ] =realengine( PI,TAU,T0,P0,M0,T,P,CP,lamb,alpha)
%Inputs necessaris: PI.c (compresor+fan) PI.f, alpha, T0, a0, P0,
%rho0,v0,M0,T.t0,P.t0,TAU.lamb,TAU.r,Tt4, CP (hot & cold), lamb (hot&cold)
%%Funció per a computar els parámetres desde l'entrada d'aire fins a la
%%sortida de la turbina d'un motor real. 
PI.cH=PI.c/PI.f;
%Eficiencies rati de pressió:
h=43e6; 
PI.d=0.96;
ETA.f=0.98; %Suposició
ETA.cH=0.98;
PI.b=0.94;
ETA.b=0.99;
ETA.tH=0.87;
ETA.tL=0.87;
PI.n=0.98;
ETA.mec=0.99;
% %Per probar el codi:
% ETA.f=1;
% PI.d=0.92;
% ETA.cH=0.86;
% PI.b=0.97;
% ETA.b=0.98;
% ETA.tH=0.92;
% ETA.tL=1;
% ETA.mec=1;

%Difusor 
T.t2=T.t0;
TAU.d=T.t2/T.t0;
P.t2=P.t0*PI.d;
%Fan
P.t25=P.t2*PI.f;
TAUi.f=pi2tau(PI.f,lamb.cold);
Ti.t25=T.t2*TAUi.f;
T.t25=T.t2+(Ti.t25-T.t2)/ETA.f;
TAU.f=T.t25/T.t2;
%Dades per després del fan al flux secundari:
P.t13=P.t25;
T.t13=T.t25;
%Compresor
P.t3=P.t25*PI.cH;
TAUi.cH=pi2tau(PI.cH,lamb.cold);
Ti.t3=T.t25*TAUi.cH;
T.t3=T.t25+(Ti.t3-T.t25)/ETA.cH;
TAU.cH=T.t3/T.t25;
%Cambra de combustió
P.t4=PI.b*P.t3;
f=(CP.hot*TAU.lamb-CP.cold*TAU.r*TAU.d*TAU.cH*TAU.f)/((ETA.b*h/T0)-CP.hot*TAU.lamb);
%Turbina alta
T.t45=T.t4-(T.t3-T.t25)*(CP.cold/(ETA.mec*(1+f)*CP.hot));
TAU.tH=T.t45/T.t4;
Ti.t45=T.t4-(T.t4-T.t45)/ETA.tH;
TAUi.tH=Ti.t45/T.t4;
PI.tH=tau2pi(TAUi.tH,lamb.hot);
P.t45=P.t4*PI.tH;
%Turbina baixa
T.t5=T.t45-(T.t25-T.t2)*(((1+alpha)*CP.cold)/((1+f)*CP.hot*ETA.mec));
TAU.tL=T.t5/T.t45;
Ti.t5=T.t45-(T.t45-T.t5)/ETA.tL;
TAUi.tL=Ti.t5/T.t45;
PI.tL=tau2pi(TAUi.tL,lamb.hot);
P.t5=P.t45*PI.tL;
end

