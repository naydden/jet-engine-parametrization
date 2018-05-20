function [ A ] = Area(M,gam,Pt,Tt,R,m)
%Cálcul área a 5. M5=1
MPF=sqrt(gam)*M*(1+((gam-1)/2)*M^2)^((-gam-1)/2*(gam-1));
A=m/MPF*sqrt(Tt*R)/Pt;
end

