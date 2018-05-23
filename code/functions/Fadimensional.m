function [ Fadim] = Fadimensional( f,M9,M19,alpha,T,P,gam,mixer,T0,M0,P0)
if mixer == true
    gam.hot=gam.mixer;
    M19=M0;
    T.s19 = T.s9;
    P.s19 = P.s9;
end

Fadim=(1+f)*M9*sqrt(T.s9/T0)-M0+alpha*(M19*sqrt(T.s19/T0)-M0)+...
    (1+f)/(M9*gam.hot)*sqrt(T.s9/T0)*(1-P0/P.s9)+alpha/(M19*gam.cold)*sqrt(T.s19/T0)*(1-P0/P.s19);

end

