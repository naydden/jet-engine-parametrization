function [ P,f_AB ] = AfterBurner( PI,P,CP,TAU,gam,ETA,h,T0,R, T,f)
%Pg 455 pdf elements of gass turbines
%After Burner parameters
gam.AB = gam.hot;
CP.AB = CP.hot;
T.t7 = T.t4;
ETA.AB = ETA.b;

%compute main parameter
R.AB = (gam.AB - 1)/gam.AB*CP.AB;
TAU.lambAB = CP.AB*T.t7/(CP.cold*T0);
%fuel fraction afterburner definded on pg 453. f_AB = m_fAB/m0
f_AB = (1+f)*(TAU.lambAB-TAU.lamb*TAU.t)/(ETA.AB*h/(CP.cold*T0)-TAU.lambAB);


end

