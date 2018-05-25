function [phi] = calcPhi(M, gamma)
phi = M^2*(1+((gamma-1)/2)*M^2)/(1+gamma*M^2)^2;
end