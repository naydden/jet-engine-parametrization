function [ P,TAU,T ] = Fan(P,PI,gam,ETA,T,TAU)
%IN: P.t2, PI.f, gam, ETA.f
%OUT:T.t13, P.t13, T.t25, P.t25, TAU.f

%Flux primari
P.t25=P.t2*PI.f;
TAUi.f=pi2tau(PI.f,gam.cold);
Ti.t25=T.t2*TAUi.f;
T.t25=T.t2+(Ti.t25-T.t2)/ETA.f;
TAU.f=T.t25/T.t2;
%Flux secundari
P.t13=P.t25;
T.t13=T.t25;
end