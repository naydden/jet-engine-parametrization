function [ T,P,M9 ] = ToveraPrimari( P0,T,PI,gam,P,TAU,mixer )
T.t9=T.t6;
P.t9=PI.n*P.t6;
P.s9=P0;
if mixer == true
    gam.hot=gam.mixer;
end
M9=sqrt((2/(gam.hot-1))*((TAU.r*TAU.d*TAU.c*TAU.b*TAU.t*TAU.n)-1));
if M9>1
    M9=1;
    P.s9=P.t9/(1+((gam.hot-1)/2)^(gam.hot/(gam.hot-1)));
end
T.s9=T.t9/(1+((gam.hot-1)/2)*M9^2);

end

