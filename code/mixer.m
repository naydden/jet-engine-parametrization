function [M6, Tt6, Pt6] = mixer(Tt5, Pt5, Tt13, Pt13, Cp5, Cp13, alpha, f)

% this function computes the mixer of a turbofan.

alpha_p = alpha/(1+f); % alpha_p = m_sec/m_prim
Cp6 = (Cp5 + alpha_p*Cp13)/(1+alpha_p);
Tt6 = (Cp5*Tt5 + alpaha_p*Cp13*Tt13)/((1+alpha_p)*Cp6);

% Say M5 = 1;
M5 = 1;
% Say gamma = 1.4
gamma = 1.4;
P5 = Pt5 / ((1+gamma/2)^(gamma/(gamma-1)));

% Since we want the flux at 6 to be as less turbulent as possible, we say
% that P5 = P13 = P6
P13 = P5;
P6 = P5;

M13 = sqrt(2/(gamma-1)*((Pt13/P13)^((gamma-1)/(gamma)))-1);

% we find M6 now by computing F6 = F5 + F13
phi_6 = R*Tt6/(gamma*((sqrt(R*Tt5/(gamma*phi(M5,gamma)))+alpha_p*sqrt(R*Tt13/(gamma*phi(M13,gamma))))/(1+alpha_p))^2);
f = @(y) phi_F(y,gamma,phi_6);
M6 = fsolve(f,0);

Pt6 = P6*(1+(gamma-1)/2*M6)^(gamma/(gamma-1));



end

