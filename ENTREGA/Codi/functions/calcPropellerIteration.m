function [Tcore, Tprop, Pow, D, n, CThrust, CPower, Jprop ] = calcPropellerIteration(m0, CP, C, v0, T0, ETA, T, rho0)
% tabla de valores passo40º
%CP     CT      J
itVal = [0.3    0.175   0.3
    0.29   0.17    0.5
    0.28   0.165   0.65
    0.27   0.16    0.8
    0.265  0.155   0.85
    0.26   0.15    0.9
    0.24   0.14    1.3
    0.23   0.13    1.4
    0.22   0.12    1.5
    0.21   0.11    1.55
    0.19   0.1     1.6
    0.175  0.09    1.65
    0.165  0.08    1.7
    0.145  0.07    1.8
    0.135  0.06    1.85
    0.12   0.05    1.9
    0.11   0.04    1.95
    0.1    0.03    2
    0.06   0.02    2.1
    0.05   0.01    2.15
    ];
[itI, itJ] = size(itVal);
%thrust

Tcore = C.cin*CP.hot*T0/v0*m0; %Tcore
Tprop = C.prop*CP.hot*T0/v0*m0*ETA.mec; %Tprop

%power
Pow = ETA.prop*ETA.mec*m0*CP.hot*(T.t45-T.t5);

%iteration parameters
finD = 4;
incD = 0.1;
finn = 170;
inc_n = 5;
D = 0.1:incD:finD; %Diametres
n = 1:inc_n:finn; %velocitats
%valors objectiu del gràfic
TOL = 0.09;
finish = false;
for diam = D
    for rev = n
        CTit = Tprop/(rho0*rev^2*diam^4);
        CPit = Pow/(rho0*rev^3*diam^5);
        Jit = v0/(rev*diam);
        for i = 1:itI
            CTref = itVal(i, 2);
            CPref = itVal(i, 1);
            Jref = itVal(i, 3);
            %check convergence
            if (abs(CTit - CTref)<TOL && abs(CPit - CPref)<TOL && abs(Jit - Jref)<TOL)
                D = diam;
                n = rev;
                CThrust = CTit;
                CPower = CPit;
                Jprop = Jit;
                finish = true;
                break;
            end
            if finish
                break;
            end
        end
        if finish
            break
        end
    end
    if finish
        break;
    end
end

end