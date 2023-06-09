%% TAKEOF-WEIGHT BUILDUP 👍😎
clc; clear all; close all; format compact;

% Crew-Weight
Wcrew = 2*(85+25); % [ kg ]  % Obtido com base em estimativa de Peso Popular

% Payload-Weight
Wpayload = 1010;   % [ kg ]  % Com base nos tres primeiros 9*(85+25) Peso Popular + Peso Bagagem 

% Altitudes de Operação
hc = 45000 * .3048;    % Teto operacional da aeronave [m] {ESTIMADO DADOS DA TABELA}
hloit = 25000 * .3048; % Altitude de Loiter [m] {ESTIMADO DADOS DA TABELA} 

% Atmosfera Padrão
[TFL450, ~, ~, ~] = atmosisa(hc);
[TFL250, ~, ~, ~] = atmosisa(hloit);

% Specific Fuel Consumption C or SFC {RAYMER TABLE 3.3 P.66} % Utilizar como base na analise de mercado dos motores utilizados nas aeronaves
SFC.FL0    =  .44/3600; % [ 1/s ]  Obtido com base nos dados da tabela
SFC.FL450  =  .5/3600; %SFC.FL0 * sqrt( TFL450 / 288.15 ); % SFC.CRUISE [ 1/s ] {AIRCRAFT PERFORMANCE SADRAEY P.150 EQ.4.21-4.30}
SFC.FL250  =  .4/3600; %SFC.FL0 * sqrt( TFL250 / 288.15 ); % SFC.LOITER [ 1/s ] {AIRCRAFT PERFORMANCE SADRAEY P.150 EQ.4.21-4.30}

% Lift-to-Drag L/D % ==================================================================================|
Cd_0 = .02; k = .0673; AR = 9;        % ======== UPDATE FROM GUDMUNDSSON                               |
LDMAX = 1 / ( 2 * sqrt( Cd_0 * k ) ); % {AIRCRAFT PERFORMANCE SADRAEY P.255 EQ.6.9}                    | 
                                      % {AIRCRAFT DESIGN SADRAEY P.226 EQ.5.21}                        |
                                      % {AIRCAFT DESIGN SADRAEY TABLE 4.5 P.131 L/D - {12-20}-{10-15}} |
                                      % {AIRCRAFT DESIGN SADRAEY P.228 TABLE 5.8 AR}                   |
% =====================================================================================================|

fprintf( '\n ============================ \n')
fprintf( '\n      A RAZÃO LD MAX É \n    LDMAX = %G \n', LDMAX ) % [ dia ]


% Cruise set-up analysis
CLB.ANG = 12; % [ deg ] Angulo de subida {!REF} {CONSULTAR MERCADO}
R.CLB   = hc / tand(CLB.ANG); % [ m ] Distancia de subida {!REF} 

    
% Loiter set-up analysis
CL.MD  = sqrt( Cd_0 / k ); % Coeficente de sustentação de minimo arrasto {!REF}
V.LOIT = ( 850 / 3.6 ) / ( 3 ^ ( 1 / 4 ) ); % Velocidade de Loiter [m/s] {!REF}   
V.MD   = V.LOIT;   % Velocidade de minimo arrasto {!REF}
E      =  45 * 60; % Tempo de espera [ s ] 
R.LOIT  = V.LOIT * E; % Distancia de Loiter [ m ] 

% Fuel-Fraction {RAYMER P.64}
Wtakeoff = .970; % W1/W0 % Sadraey define como 0.980 {SADRAEY TABLE 4.3 P.129}
Wclimb   = .985; % W2/W1 % Sadraey define como 0.970 {SADRAEY TABLE 4.3 P.129}
Wdescent = .990; % W5/W4 % Sadraey define como 0.990 {SADRAEY TABLE 4.3 P.129}
Wlanding = .995; % W6/W5 % Sadraey define como 0.997 {SADRAEY TABLE 4.3 P.129}

% Cruise-Fuel Fraction {RAYMER EQ.3.6} ================================================================== |
R      = 3900*1E3 - R.CLB - R.LOIT; %(3750*1.40)*1E3 - R.CLB - R.LOIT ; % Range [ m ] {BASE CONCORRENTES}                        |
V.CRU  = 850 / 3.6; % Velocidade de Cruseiro [ m/s ]                                                      |
Wcruise = exp( ( -R * SFC.FL450 ) / ( V.CRU * (.866 * LDMAX) ) ); % W3/W2 {3.4.5 P.71}                    |
%                                                                                                         |
fprintf( '\n ============================ \n ') %                                                         |   
fprintf( '\n   A DISTANCIA PERCORRIDA EM \n') %                                                           |
fprintf( '     CRUZEIRO DA AERONAVE É \n          R = %G \n', R * 1E-3 ) %                                |
%                                                                                                         |
fprintf( '\n   AUMENTO DE %G km COM  \n', R * 1E-3 - 3750) %                                          |
fprintf( '     RELAÇÃO A REFERÊNCIA \n') %                                                                        |
% ======================================================================================================= |

% Loiter-Fuel Fraction {RAYMER EQ.3.6} ================================================================== |
Wloiter = exp( ( -E * SFC.FL250 ) / ( LDMAX ) ); % W4/W3 {SADRAEY P.133 EQ.4.20} |
% ======================================================================================================= |

% Overall fuel weight ratio
WxW0 = Wtakeoff * Wclimb * Wcruise * Wloiter *  Wlanding * Wdescent; % W6/W0 {SADRAEY P.128 EQ.4.10}

% fuel weight ratio
W.f0 = 1.06 * ( 1 - WxW0 ); % 5 a mais de combustível {RAYMER P.71} {SADRAEY P.129 EQ.4.11}

% Takeof-Weight Guess
W0G = linspace( 2E3, 2E5, 1E6 );

% Empty-Weight {RAYMER P.59} {SADRAEY P.136 EQ.4.26}
W.e0 =1.4 * W0G .^ ( -.1 ); % {RAYMER P.61 TABLE.3.1}
% W.e0 = -7.754*1E-8 * W0G + .576; % {SADRAEY P.136/138 EQ.4.26 TABLE 4.8} 

% Takeof-Weight
W0  = ( Wcrew + Wpayload )./ ( 1 - W.f0 - W.e0 ); % Descobre MTOW {SADRAEY P.123 EQ.4.5}

% DIFERENCE
WDif = W0-W0G;
[~,idx] = min(abs(W0-W0G));
WDif = WDif(1:idx);


%% GRÁFICOS

figure
yyaxis left
plot(W0G(1:length(WDif)),W0(1:length(WDif)), 'LineWidth', 2)
ylabel('Takeoff-Weight (W0) [kg]')
yyaxis right
plot(W0G(1:length(WDif)),WDif, '--','LineWidth', 1)
xline(W0G(length(WDif)),'-.','MAX WEIGHT')
xlabel('Takeoff-Weight Guess (W0) [kg]')
xlim ( [ W0G(1) W0G(length(WDif))+.5E3 ] )
title('Takeoff-Weight Sizing')
legend('Takeoff-Weight',...
       'Difference',...
       'FontSize', 12,'FontName','Times New Roman','Location','best');
grid on; grid minor;
set(gcf, 'Color', 'w');
set(gca,'GridLineStyle', '-');
set(gcf,'paperPositionMode','auto')  


% DEFINE PESO
W0 = round(W0G(idx));

fprintf( '\n ============================ \n')
fprintf( '\n     O peso da aeronave é \n') % 
fprintf( '            W = %G \n', W0 ) % 
fprintf('\n ======== FIM %s ============= \n',char( double( '👍😎' ) ) )

% =======================================================================
