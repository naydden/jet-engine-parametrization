function [ PI ] = tau2pi( TAU,lamb )
%tau2pi From tau to pi values
PI = TAU^(lamb/(lamb-1));


end

