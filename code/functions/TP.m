function [C, ETA] = TP(P,PI,gam,ETA,T,TAU,T0,M0,f,M9,P0, R, h, CP )
%pg507 ELEMENTS OF GAS TURBINES
ETA.g = 0.99;
ETA.mL= 0.99;

C.cin = (gam.cold-1)*M0*((1+f)*M9-M0+(1+f)*R.hot/R.cold*((T.s9/T0)/M9)*(1-P0/P.s9)/gam.cold);
C.prop = ETA.prop*ETA.g*ETA.mL*(1+f)*TAU.lamb*TAU.tH*(1-TAU.tL);
C.tot = C.cin + C.prop;

%propulsive efficiency
ETA.P = C.tot/(C.prop/ETA.prop+((gam.cold-1)/2)*((1+f)*M9^2-M0^2));
ETA.T = C.tot/((f*h)/CP.cold*T0);


end