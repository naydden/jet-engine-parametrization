function [ T,P,M9,PI,TAU ] = ToveraPrimari( P0,T,PI,gam,P,TAU,isMixer,isTurboProp,ETA,isAfterBurner )
T.t9=T.t6;
% if isAfterBurner
%     T.t9=T.t7;
%     gam.hot=gam.AB;
% end
P.t9=PI.n*P.t6;
P.s9=P0;
if isMixer
    gam.hot=gam.mixer;
end
if isTurboProp
    TAU.f=1;
end
M9=sqrt((2/(gam.hot-1))*((TAU.r*TAU.d*TAU.cH*TAU.f*TAU.b*TAU.t*TAU.n)-1));
if ~isreal(M9)
    M9=0.1; %Fixar M9 a un valor baix
    %Considerar tovera adaptada. Es fa el procés invers per saber la mínima
    %pressió de sortida a la turbina. 
    P.t9=P.s9*(1+(gam.hot-1)/2*M9^2)^gam.hot/(gam.hot-1);
    P.t6=P.t9/PI.n;
    T.t6=T.t9;
    P.t5=P.t6;
    PI.tL=P.t5/P.t45;
    TAUi.tL=pi2tau(PI.tL,gam.hot);
    Ti.t5=TAUi.tL*T.t45;
    T.t5=T.t45-ETA.tL*(T.t45-Ti.t5);
    TAU.tL=T.t5/T.t45;
    TAU.t=TAU.tL*TAU.tH;
end
if M9>1
    M9=1;
    P.s9=P.t9/(1+((gam.hot-1)/2))^((gam.hot-1)/gam.hot);
end
T.s9=T.t9/(1+((gam.hot-1)/2)*M9^2);

end

