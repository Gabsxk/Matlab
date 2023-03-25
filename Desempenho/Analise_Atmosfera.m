%% Deempenho - 07/03/2023
clc; close all; clear all; format compact
%% Segundo Exercício
% Considere uma aeronave voando em um nível FL300, em uma atmosfera
% padrão. Determine as propriedades da atmosfera durante o voo.

% FL300 = 30 000 pés -> Troposfera

% Inicialização de Variáveis
H  = 30000*0.3048;         % Altitude [ m ]
Hv  = linspace(0,H*2,100); % Altitude variando [m]
p0 = 101325;               % Pressão [ N/m^2 ] 
L0 = -0.0065;              % Taxa de alteração da temperatura em função de H (coef. ang.) [ K/ft ]
T0 = 288.15;               % Temperatura [ K ]
g0 = 9.80665;              % Aceleração da gravidade [ m/s^2 ]
R  = 287.05;               % Constante universal dos gases [ Nm/kgK ]

%% Utilizando atmosisa
[T, a, P, rho] = atmosisa(H);
    fprintf("Temperatura = %e [K] \n", T)
    fprintf("Velocidade do som = %e [m/s] \n", a)
    fprintf("Pressão = %e [Pa]\n", P)
    fprintf("Densidade = %e [kg/m^3] \n", rho)

[Tv, av, Pv, rhov] = atmosisa(Hv);
    figure
    plot(rhov,Hv);title('Densidade-Altitude')
    xlabel('Densidade [kg/m^3]'); ylabel('Altitude [m]')
    legend('Densidade',...
    'Location','NorthWest','FontSize', 12,'FontName','Times New Roman','Location','best');
    set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
    set(gcf,'paperPositionMode','auto');
    
%% Terceiro Exercício
fprintf("\n Terceiro Exercício")
% As propriedades da atmosfera local do aeroporto de Uberlândia são:
%   • Altitude real: 3094 ft
%   • Pressão atmosférica: 95 kPa
%   • Temperatura ambiente: 27 °C
% A partir destas condições, determinar se a atmosfera encontra-se em sua
% forma padrão ISA. Caso não, determinar a altitude pressão e altitude
% densidade para um voo neste aeroporto sob tais condições.

% Verificar se para essa altitude a pressão apresentada corresponde a
% pressão na atmosfera padrão ISA 

H3 = 3094 * .3048; % Altitude [ m ]
T = 27 + 273.15;  % Temperatura [ K ]
p = 95000;        % Pressão [ Pa ] 

[Tisa, ~, P, rhoisa] = atmosisa(H3);
fprintf("\n A temperatura deveria ser 27°C e o obtido foi %G°C", Tisa-273.15)

% Pela temperatura observa-se que o aeroporto esta em uma condição ISA+18°C
% pois a temperatura padrão ISA seria 9°C e a diferença foi de (T)27 - (T_ISA)9 = 18
% °C. Assim analiza-se a atmofera para a temperatura de 18°C.
[~, index] = min(abs(T-Tisa));

rho3 = rho(index);
fprintf("\n A densidade para ISA+18 %G", rho3)



