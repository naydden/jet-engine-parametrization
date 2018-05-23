% Testing mixer.m with activity 2.5
%%  DATA
f = 0.013;
m0 = 26.2;

T0 = 288.15; %at SL
P0 = 101325; %at SL

M5 = 0.41;
Tt5 = T0*4.1;
Pt5 = P0*4.2;

M13 = 0.39;
Tt13 = T0*1.51;

m5_m6 = 0.72;
m13_m6 = 0.28;
alpha_p = (m13_m6*m0)/(m5_m6*m0);

R5 = 287;
R13 = 287;
gamma5 = 1.3;
gamma13 = 1.4;

%% CALCULATIONS
Cp13 = R13*gamma13/(gamma13-1);
Cp5 = R5*gamma5/(gamma5-1);

[M6, Tt6, Pt6] = mixer(Tt5, Pt5, Cp5, gamma5, Tt13, M13, Cp13, gamma13, alpha_p,M5)
