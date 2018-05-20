function [ P9,T9,P19,T19,Fadim,M9,M19 ] = Nozzle( T,P,PI,ETA,P0,T0,alpha,gam,f,M0)
T.t9=T.t6;
P.t9=PI.n*P.t6;
%Flux primari
P9=P0;
M9=sqrt((2/(gam.hot-1))*((PI.r*PI.d*PI.c*PI.b*PI.t*PI.n)^((gam.hot-1)/gam.hot)-1));
if M9>1
    M9=1;
    P9=P.t9/(1+((gam.hot-1)/2)^(gam.hot/(gam.hot-1)));
end
T9=T.t9/(1+((gam.hot-1)/2));
%Flux secundari
T.t19=T.t16;
P.t19=PI.n*P.t16;
P19=P0;
M19=sqrt((2/(gam.cold-1))*((PI.r*PI.d*PI.f*PI.n)^((gam.cold-1)/gam.cold)-1));
if M19>1
    M19=1;
    P19=P.t19/(1+((gam.cold-1)/2)^(gam.cold/(gam.cold-1)));
end
T19=T.t19/(1+((gam.cold-1)/2));
Fadim=(1+f)*M9*sqrt(T9/T0)-M0+alpha*(M19*sqrt(T19/T0)-M0)+(1+f)/(M9*gam.hot)*sqrt(T9/T0)*(1-P0/P9)+alpha/(M19*gam.cold)*sqrt(T19/T0)*(1-P0/P19);
end

