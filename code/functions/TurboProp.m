function [C] = TurboProp(P,PI,gam,ETA,T,TAU,T0, M0 )
%TurboProp calculation segons EXERCICI 33

C.cin = (gam.cold-1)*M0*(sqrt(2/(gam.cold-1)*(TAU.tL*(TAU.lamb+TAU.r-TAU.r*TAU.cH)-TAU.lamb/(TAU.r*TAU.cH)))-M0);
C.prop = ETA.prop*ETA.mec*(1-TAU.tL)*(TAU.lamb+TAU.r-TAU.r*TAU.cH);
C.tot = C.cin + C.prop;


end