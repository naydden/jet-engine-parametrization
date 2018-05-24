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
isMixer = false; % si esta en true col·loca el mixer
isAftBurner = false; % si esta en true calcula l'after burner
isTurboProp = false; % si esta en true col·loca un turboprop
isTP = false;
%% PRE-PROCESSING - Find design optimum parameters PI.f, PI.c, alpha
[ PI, alpha ] = opt_parameters( M0, a0, gam, gc, PI, TAU);

fprintf('The optimum values are: \n pi_f = %.2f\n pi_c = %.2f\n alpha = %.2f\n',...
    PI.f, PI.c,alpha);
%% PROCESSING - Main code
%Calcul de les etapes del jet:
[T,P,TAU] = Difusor( T,P,TAU,PI );
if isTurboProp %optimitzaci� de la turbina
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
    Tcore = C.cin*CP.hot*T0/v0; %Tcore/m0
    Tprop = C.prop*CP.hot*T0/v0; %Tprop/m0
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