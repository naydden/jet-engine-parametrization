clc; clear; close all;
%% TURBOFAN PARAMETRIC CALCULATION
%{ 20/05/2018: Pol Fontanes, Eva Maria Urbano Gonzalez i Boyan Naydenov%}
%% INPUT - carregar parametres generals
%Add to path functions folder
addpath(genpath('./functions'))
% Define input data files
NAME_INPUT_DATA = 'DATA' ; %Fitxer amb dades generals donades
% Load data files
eval(NAME_INPUT_DATA); %Carregar el codi
%% Selector seccions - INPUT
isMixer = false; % si esta en true colÂ·loca el mixer
isAftBurner = false; % si esta en true calcula l'after burner
isTurboProp = true; % si esta en true colÂ·loca un turboprop
isTP = false;
%% PRE-PROCESSING - Find design optimum parameters PI.f, PI.c, alpha
[ PI, alpha ] = opt_parameters( M0, a0, gam, gc, PI, TAU );
fprintf('The optimum values are: \n pi_f = %.2f\n pi_c = %.2f\n alpha = %.2f\n',...
    PI.f, PI.c,alpha);


%% PROCESSING - Main code
%Calcul de les etapes del jet:
[T,P,TAU] = Difusor( T,P,TAU,PI );
if isTurboProp %optimitzaciï¿½ de la turbina
    %exercici 33 - es va seguint pas a pas
    [PI,P,TAU,T] = Compressor( PI,P,gam,T,TAU,ETA,isTurboProp);
    [P,f] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0,isTurboProp );
    [T,TAU,PI,P] = TurbinaAlta( T,CP,ETA,f,gam,TAU,P,PI);
    [T,TAU,PI,P ] = TurbinaBaixa(T,alpha,CP,f,ETA,gam,P,PI,TAU,isTurboProp);
    J = 0.76;
    D = 2.6;
    n = calcPropeller(ETA, J, v0, D);
elseif isTP %afegir una propeller
    [PI,P,TAU,T] = Compressor( PI,P,gam,T,TAU,ETA,1);
    [P,f] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0,1 );
    [T,TAU,PI,P] = TurbinaAlta( T,CP,ETA,f,gam,TAU,P,PI);
    TAU.tL = 0.8751;
    TAU.t=TAU.tH*TAU.tL;
    [C] = TurboProp(P,PI,gam,ETA,T,TAU, T0, M0 );
    %         [C, ETA] = TP(P,PI,gam,ETA,T,TAU,T0,M0,f,M9,P0,R,h, CP, isTP );
    Fadim_TP = C.tot*CP.cold*T0/(v0*a0);
    Tcore = C.cin*CP.hot*T0/v0; %Tcore/m0
    Tprop = C.prop*CP.cold*T0/v0; %Tprop/m0
    J = 0.76;
    D = 2.6;
    n = calcPropeller(ETA, J, v0, D);
else %cas estandard turbofan
    [P,TAU,T] = Fan( P,PI,gam,ETA,T,TAU );
    [PI,P,TAU,T] = Compressor( PI,P,gam,T,TAU,ETA,isTurboProp);
    [P,f] = CambraCombustio( PI,P,CP,TAU,gam,ETA,h,T0,isTurboProp );
    [T,TAU,PI,P] = TurbinaAlta( T,CP,ETA,f,gam,TAU,P,PI);
    [T,TAU,PI,P ] = TurbinaBaixa(T,alpha,CP,f,ETA,gam,P,PI,TAU,isTurboProp);
end
%%
%Punts d'entrada a mixer: 1.3 i 5. Punt a la sortida del mixer: 6
if isMixer
    %COMPUTE THE MIXER
    [ T,P,M6,gam,CP,R] = mixer(T,P,CP,gam,R,alpha,f,M0);
    %[ T, P, M6A, gam, CP] = mixer2(T,P, CP,gam,R, alpha,f, T0, PI, TAU);
    TAU.b=tau2pi(PI.b,gam.mixer);
    TAU.n=tau2pi(PI.n,gam.mixer);
    [ T,P,M9,PI,TAU ] = ToveraPrimari( P0,T,PI,gam,P,TAU,isMixer,isTurboProp,ETA);
    Fadim = Fadimensional( f,M9,M9,alpha,T,P,gam,isMixer,T0,M0,P0,isAftBurner);
    
elseif isTurboProp %Tovera si tenim turboprop i no tenim mixer
    P.t6=P.t5;
    T.t6=T.t5;
    TAU.b=tau2pi(PI.b,gam.hot);
    TAU.n=tau2pi(PI.n,gam.hot);
    [ T,P,M9,PI,TAU ] = ToveraPrimari( P0,T,PI,gam,P,TAU,isMixer,isTurboProp,ETA,isAftBurner);
    [C] = TurboProp(P,PI,gam,ETA,T,TAU, T0, M0 );
%     [C, ETA] = TP(P,PI,gam,ETA,T,TAU,T0,M0,f,M9,P0,R,h,CP, isTP );
    Fadim_TP = C.tot*CP.cold*T0/(v0*a0);
        
    %proppeller iteration
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
    [ m0,mf,msec ] = Fluxosmasics( f,Fadim_TP,F,a0,0);
    Tcore = C.cin*CP.hot*T0/v0*m0; %Tcore
    Tprop = C.prop*CP.hot*T0/v0*m0*ETA.mec; %Tprop
    
        %power
    Pow = ETA.prop*ETA.mec*m0*CP.hot*(T.t45-T.t5);
        
    %iteration parameters
    finD = 2; 
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
    
elseif isTP

else %Tovera sense turboporp ni mixer
    P.t6=P.t5;
    T.t6=T.t5;
    TAU.b=tau2pi(PI.b,gam.hot);
    TAU.n=tau2pi(PI.n,gam.hot);
    P.t16=P.t13;
    T.t16=T.t13;
    [ T,P,M9,PI,TAU ] = ToveraPrimari( P0,T,PI,gam,P,TAU,isMixer,isTurboProp,ETA,isAftBurner);
    [ T,P,M19] = Toverasecundari( T,PI,P,P0,gam,TAU );
    Fadim = Fadimensional( f,M9,M19,alpha,T,P,gam,isMixer,T0,M0,P0,isAftBurner);
    if isTP %afegir una propeller 
%         [C] = TurboProp(P,PI,gam,ETA,T,TAU, T0, M0 );
% %         [C, ETA] = TP(P,PI,gam,ETA,T,TAU,T0,M0,f,M9,P0,R,h, CP, isTP );
%         Fadim_TP = C.tot*CP.cold*T0/(v0*a0);
%         Tcore = C.cin*CP.hot*T0/v0; %Tcore/m0
%         Tprop = C.prop*CP.hot*T0/v0; %Tprop/m0
%         J = 0.76;
%         D = 2.6;
%         n = calcPropeller(ETA, J, v0, D);
    end
end
%Afegir afterburner
if isAftBurner
    [ P,f_AB, fab, Fadim_prim_AB, T] = AfterBurner( PI,P,CP,TAU,gam,ETA,h,T0,R, T,f, M9, M0, P0);
    f = f + f_AB;
    [ T,P,M9,PI,TAU ] = ToveraAFT( P0,T,PI,gam,P,TAU,ETA);
    %Empenta adimensional total
    Fadim_AB = Fadimensional( f,M9,M19,alpha,T,P,gam,isMixer,T0,M0,P0,isAftBurner);
    [ m0af,mfaf,msecaf ] = Fluxosmasics( f,Fadim_AB,F,a0,alpha);
    
end

if ~isTurboProp && ~isTP
    [ m0,mf,msec ] = Fluxosmasics( f,Fadim,F,a0,alpha);
    
    %Calcul Arees
    %Fluxos masics:
    if isAftBurner
        m5=m0af+mfaf;
        msortida=m0af+mfaf+msecaf;
        mentrada=m0af+msecaf;
    else
        m5=m0+mf;
        msortida=m0+mf+msec;
        mentrada=m0+msec;
    end
    
    A.e0 = Area(M0,gam.cold,P.t0,T.t0,R.cold,mentrada);
    if isMixer==true
        A.e9 = Area(M9,gam.mixer,P.t9,T.t9,R.mixer,msortida);
    else
        A.e9 = Area(M9,gam.hot,P.t9,T.t9,R.hot,m5);
        A.e19 = Area(M19,gam.cold,P.t19,T.t19,R.cold,msec);
    end
    
end