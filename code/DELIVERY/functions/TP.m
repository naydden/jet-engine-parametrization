function [C, ETA] = TP(P,PI,gam,ETA,T,TAU,T0,M0,f,M9,P0, R, h, CP, isTP )
%pg507 ELEMENTS OF GAS TURBINES
if isTP
    P.t9_P0 = PI.r*PI.d*PI.cH*PI.b*PI.t*PI.n;
else
    P.t9_P0 = P.t9/P0;
end
if P.t9_P0 > ((gam.hot+1)/2)^(gam.hot/(gam.hot-1))
    M9 = 1;
    P.t9_9 = ((gam.hot+1)/2)^(gam.hot/(gam.hot-1));
    P.s0_s9 = P.t9_9/P.t9_P0;
else
    M9 = sqrt(2/(gam.hot-1)*((P.t9_P0)^((gam.hot-1)/gam.hot)-1));
    P.t9_9 = P.t9_P0;
    P.s0_s9 = 1;
end
V.v9_a0 = sqrt(2*TAU.lamb*TAU.tH*TAU.tL/(gam.cold-1)*(1-(P.t9/P.s9)^-((gam.hot-1)/gam.hot)));
%Ho calculo jo T9/T0 
%T.s9_T0 = 1/(1+(gam.hot-1)/2*M9^2)*TAU.n*TAU.t*TAU.b;

C.cin = (gam.cold-1)*M0*((1+f)*V.v9_a0-M0+(1+f)*R.hot/R.cold*((T.s9/T0)/V.v9_a0)*(1-P.s0_s9)/gam.cold);
C.prop = ETA.prop*ETA.g*ETA.mL*(1+f)*TAU.lamb*TAU.tH*(1-TAU.tL);
C.tot = C.cin + C.prop;

%propulsive efficiency
ETA.P = C.tot/(C.prop/ETA.prop+((gam.cold-1)/2)*((1+f)*M9^2-M0^2));
ETA.T = C.tot/((f*h)/CP.cold*T0);



end