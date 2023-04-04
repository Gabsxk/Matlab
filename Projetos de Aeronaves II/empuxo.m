%% Calculo do empuxo 
%potencia disponivel - potencia requerida
format compact
clc; clear; close all;
T_max             = 2*[23219.7084	   15470.90916	     31137.54	   16013.592	     16458.414	    20808.77316];
%                   *   PW306B	     PW535E	           LF507	     FJ44-4	         TFE-731-2	       pw305A
h                 = 45000*0.3048;
[Temp, a, P, rho] = atmosisa(h);
Tfl               = T_max*((rho/1.225))^1;
V                 = [1:1:500];
ew                = interp1([0 30],[.7831 .5071],[28],'linear'); 
ARw               = 9;                   
k2                = 1/(pi*ARw*ew);       

Sw                = 25.6858;
w                 = 6642*9.81; 
Cfe               = 0.003;               
 
Swet_Sref         = 6;                   

Cd0               = Cfe*Swet_Sref;    
 

Pr                = 0.5*rho*(V.^3)*Sw*Cd0 + (2*k2*w^2)./(rho*V.*Sw);
for ii=1:length(Tfl)
    Pd                = Tfl(ii).*V;
    delta_P = Pd - Pr;
    figure
    hold on
    plot(V,Pr)
    plot(V,Pd)
    plot(V,delta_P)
    [~,index] = max(delta_P);
    plot(V(index),delta_P(index),'o')
    ylim([0 3E6])
    xlim([50 400])
    ROC               = delta_P(index)/w; %rage of climb
    mask_results_tfl(ii) = ROC>=2.54;
end
mask_results_tfl
%T = (2.54*w + Pr)/V