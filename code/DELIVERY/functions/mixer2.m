function [ T, P, M6A, gam, CP] = mixer2(T,P, CP,gam,R, alpha,f, T0, PI, TAU)
% this function computes the mixer of a turbofan. PG 490
alpha_p = alpha/(1+f);
R.mixer = (R.hot + alpha_p*R.cold)/(1+alpha_p);
CP.mixer = (CP.hot + alpha_p*CP.cold)/(1+alpha_p);
gam.mixer = CP.mixer/(CP.mixer - R.mixer);

M6 = 1;

%P.t16_t6 = PI.f/(PI.c*PI.b*PI.t);
P.t16_t6 = 1;
M16 = sqrt(2/(gam.cold-1)*((P.t16_t6*(1+(gam.hot-1)/2*M6^2)^(gam.hot/(gam.hot-1)))^((gam.cold-1)/gam.cold)...
    -1));

T.t16_t6 = T0*TAU.r*TAU.f/(T.t4*TAU.t);

TAU.M = (CP.hot*(1+alpha_p*(CP.cold/CP.hot)*(T.t16_t6)))/(CP.mixer*(1+alpha_p));

phi6 = calcPhi(M6, gam.hot);
phi16 = calcPhi(M16, gam.cold);

PHI = ((1+alpha_p)/(1/sqrt(phi6)+alpha_p*sqrt(R.cold*gam.hot*T.t16_t6/(...
    R.hot*gam.cold*phi16))))^2*R.mixer*gam.hot/(R.hot*gam.mixer)*TAU.M;

M6A = sqrt(2*PHI/(1-2*gam.mixer*PHI+sqrt(1-2*(gam.mixer+1)*PHI)));

% S'ha de revisar
T.t6 = (CP.hot*T.t5 + alpha_p*CP.cold*T.t13)/((1+alpha_p)*CP.mixer);

end

