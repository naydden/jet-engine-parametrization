function [ T,P,TAU ] = Difusor( T,P,TAU,PI )
%IN: T.t0, P.t0, PI.d; 
%OUT: T.t2, TAU.d, P.t2
T.t2=T.t0;
TAU.d=T.t2/T.t0;
P.t2=P.t0*PI.d;
end

