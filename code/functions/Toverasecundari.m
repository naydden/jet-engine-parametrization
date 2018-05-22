function [ T,P,M19] = Toverasecundari( T,PI,P,P0,gam,TAU )
T.t19=T.t16;
P.t19=PI.n*P.t16;
P.s19=P0;
M19=sqrt((2/(gam.hot-1))*((TAU.r*TAU.d*TAU.f*TAU.n)-1));
if M19>1
    M19=1;
    P.s19=P.t19/(1+((gam.cold-1)/2))^((gam.cold-1)/gam.cold);
    
end
T.s19=T.t19/(1+((gam.cold-1)/2)*M19^2);
end

