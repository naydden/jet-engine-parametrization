function [ T, P, M6, gam, CP,R] = mixer(T,P, CP,gam,R, alpha,f,M0)
% this function computes the mixer of a turbofan.
alpha_p = alpha/(1+f);
R.mixer = (R.hot + alpha_p*R.cold)/(1+alpha_p);
CP.mixer = (CP.hot + alpha_p*CP.cold)/(1+alpha_p);
gam.mixer = CP.mixer/(CP.mixer - R.mixer);

T.t6 = (CP.hot*T.t5 + alpha_p*CP.cold*T.t13)/((1+alpha_p)*CP.mixer);

% Say M5 = 1;
M5 = 1;
P.s5 = P.t5 / ((1+((gam.hot-1)*0.5*M5^2))^(gam.hot/(gam.hot-1)));

% Since we want the flux at 6 to be as less turbulent as possible, we say
% that P5 = P13 = P6
P.s13 = P.t2 / ((1+((gam.cold-1)*0.5*M0^2))^(gam.cold/(gam.cold-1))); %P2 = P13
P.s6 = P.s5;

M13 = sqrt(2/(gam.cold-1)*((P.t13/P.s13)^((gam.cold-1)/(gam.cold))-1));

% we find M6 now by computing F6 = F5 + F13
phi_6 = R.mixer*T.t6/(gam.mixer*((sqrt(R.hot*T.t5/(gam.hot*phi(M5,gam.hot)))+alpha_p*sqrt(R.cold*T.t13/(gam.cold*phi(M13,gam.cold))))/(1+alpha_p))^2);
% if phi_6 > 0.208 mixer ahogado y se supone M6 = 1;
if phi_6 > 0.208
    M6 = 1;
else
    M6s = sym('M6s');
    assume(M6s <= 1 & M6s > 0);
    eqn = phi_6 == ((M6s*sqrt(1+(gam.mixer-1)*0.5*M6s^2))/(1+gam.mixer*M6s^2))^2;
    M6 = double(solve(eqn,M6s));
end
P.t6 = P.s6*(1+(gam.mixer-1)/2*M6)^(gam.mixer/(gam.mixer-1));

end

