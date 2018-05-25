function [ T,P,M9 ] = ToveraTP( P0,T,PI,gam,P,TAU )
%Nozzle
TAU.n = 1; PI.n = 1;
P.t5 = PI.tL*P.t45;
P.t9 = P.t5;

%Exhaust
P9 = P0;
M9 = sqrt(2/(gam.hot-1)*((P.t9/P9)^((gam.hot-1)/gam.hot)-1));
%Alert through workspace.
if ~isreal(M9)
    warning('M9 imaginary; Pt9/P9 = %.2d ', P.t9/P0);
end
end

