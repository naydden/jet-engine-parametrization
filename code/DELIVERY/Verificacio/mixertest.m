function [M6, Tt6, Pt6] = mixertest(Tt5, Pt5, Cp5, gamma5, Tt13, M13, Cp13, gamma13, alpha_p, M5)
R = 287;
% this function computes the mixer of a turbofan.

% alpha_p = alpha/(1+f); % alpha_p = m_sec/m_prim
Cp6 = (Cp5 + alpha_p*Cp13)/(1+alpha_p);
Cv6 = Cp6 - R;
gamma6 = Cp6/Cv6;
Tt6 = (Cp5*Tt5 + alpha_p*Cp13*Tt13)/((1+alpha_p)*Cp6);

% Say M5 = 1;
% M5 = 1;
% Say gamma = 1.4
P5 = Pt5 / ((1+((gamma5-1)*M5^2)/2)^(gamma5/(gamma5-1)));

% Since we want the flux at 6 to be as less turbulent as possible, we say
% that P5 = P13 = P6
P13 = P5;
P6 = P5;

% M13 = sqrt(2/(gamma-1)*((Pt13/P13)^((gamma-1)/(gamma)))-1);

% we find M6 now by computing F6 = F5 + F13
phi_6 = R*Tt6/(gamma6*((sqrt(R*Tt5/(gamma5*phi(M5,gamma5)))+alpha_p*sqrt(R*Tt13/(gamma13*phi(M13,gamma13))))/(1+alpha_p))^2)

M6s = sym('M6s');
assume(M6s < 1 & M6s > 0);
eqn = phi_6 == ((M6s*sqrt(1+(gamma6-1)*0.5*M6s^2))/(1+gamma6*M6s^2))^2;
M6 = double(solve(eqn,M6s));

Pt6 = P6*(1+(gamma6-1)/2*M6)^(gamma6/(gamma6-1));



end