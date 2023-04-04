clc; clear all; close all
%% PROJETO DE AERONAVES II
%% ESTIMATIVA DA POLAR DE ARRASTO
% PAR�METROS
% bw              = ;                    % Envergadura [m]
Sw                = 20.6858;             % Area da Asa [m^2]
% Swet            = ;                    % Area da Molhada da Asa [m^2]
w                 = 6642*9.81;           % MTOW [kg] - MODIFICAAAAR
ARw               = 9;                   % Alongamento da Asa (bw^2)/Sw

Cfe               = 0.003;               % Coeficiente de fric��o (aeronave de transporte de passageiros a jato)
 
Swet_Sref         = 6;                   % Raz�o de �rea molhada Swet/Sw (aeronave de transporte de passageiros a jato = 6)

Cd0               = Cfe*Swet_Sref;       % Coeficiente de Arrasto Parasita
Swet              = 0.0317*w^(0.7530);
Cd0_teste         = Cfe*(Swet/Sw);

ew                 = interp1([0 30],[.7831 .5071],[28],'linear');           % Interpola��o para melhor obten��o de eficiencia de Evergadura
%ew                = 1/(1+VALORX);        % Constante de Oswald

k2                 = 1/(pi*ARw*ew);       % Constante do Coeficiente de Arrasto Induzido

%% FLIGHT CONDITIONS
h               = 45000*0.3048;          % Altitude
[T, a, P, rho]  = atmosisa(h);
v               = 90:1:236.1;            % Velocidade [m/s]
Cl              = [0:0.001:3];           % Coeficiente de sustenta��o

%% POLAR DE ARRASTO
%Cd = Cd0 (Arrasto Parasita) + Cdi (Arrasto Induzido)
Cd              = Cd0 + k2 * Cl.^2;
plot(Cl,Cd,'linewidth', 2)
grid minor
xlabel('C_L')
ylabel('C_D')



