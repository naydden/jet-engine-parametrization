function F = phi_F(M,gamma,phi)
F = phi - ((M*sqrt(1+(gamma-1)*0.5*M^2))/(1+gamma*M^2))^2;
end