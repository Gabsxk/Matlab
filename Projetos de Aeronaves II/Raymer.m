%% TAKEOF-WEIGHT BUILDUP 👍😎
clc; clear all; close all; format compact;

% Crew-Weight
Wcrew = 2*(85+25); % [ kg ]  % Obtido com base em estimativa de Peso Popular

% Payload-Weight
Wpayload = 1000;   % [ kg ]  % Com base nos tres primeiros 9*(85+25) Peso Popular + Peso Bagagem 

% Altitudes de Operação
hc = 45000 * .3048;    % Teto operacional da aeronave [m] {ESTIMADO DADOS DA TABELA}
hloit = 25000 * .3048; % Altitude de Loiter [m] {ESTIMADO DADOS DA TABELA} 

% Atmosfera Padrão
[TFL450, ~, ~, ~] = atmosisa(hc);
[TFL250, ~, ~, ~] = atmosisa(hloit);

% Specific Fuel Consumption C or SFC {RAYMER TABLE 3.3 P.66} % Utilizar como base na analise de mercado dos motores utilizados nas aeronaves
SFC.FL0    =  .44/3600;                          % [ 1/s ]  Obtido com base nos dados da tabela
SFC.FL450  =  SFC.FL0 * sqrt( TFL450 / 288.15 ); % SFC.CRUISE [ 1/s ] {AIRCRAFT PERFORMANCE SADRAEY P.150 EQ.4.21-4.30}
SFC.FL250  =  SFC.FL0 * sqrt( TFL250 / 288.15 ); % SFC.LOITER [ 1/s ] {AIRCRAFT PERFORMANCE SADRAEY P.150 EQ.4.21-4.30}

% Lift-to-Drag L/D
Cd_0 = .02; k = .0673; AR = 9;        % {CÓDIGO GUDMUNDSOSN}
Cd_0 = .02; k = .0551; AR = 11;       % {CÓDIGO GUDMUNDSOSN} % OSWALD NAO MUDA {AR 9} LD AUMENTA {PESO DIMINUI}
LDMAX = 1 / ( 2 * sqrt( Cd_0 * k ) ); % {AIRCRAFT PERFORMANCE SADRAEY P.255 EQ.6.9} 
                                      % {AIRCRAFT DESIGN SADRAEY P.226 EQ.5.21}
                                      % {AIRCAFT DESIGN SADRAEY TABLE 4.5 P.131 L/D - {12-20}-{10-15}}
                                      % {AIRCRAFT DESIGN SADRAEY P.228 TABLE 5.8 AR}

fprintf( '\n ============================ \n')
fprintf( '\n      A RAZÃO LD MAX É \n    LDMAX = %E \n', LDMAX ) % [ dia ]


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
Wtakeoff = .970; % W1/W0
Wclimb   = .985; % W2/W1
Wlanding = .995; % W5/W4
Wdescent = .990; % W6/W5

% Cruise-Fuel Fraction {RAYMER EQ.3.6} ===========================================|
R      = (3750*1.95)*1E3 - R.CLB - R.LOIT ;      % Range [ m ] {BASE CONCORRENTES}                             |
V.CRU  = 850 / 3.6; % Velocidade de Cruseiro [ m/s ]                              |
Wcruise = exp( ( -R * SFC.FL450 ) / ( V.CRU * (.866 * LDMAX) ) ); % W3/W2 {3.4.5 P.71}|
% ================================================================================|

fprintf( '\n ============================ \n ')
fprintf( '\n   A DISTANCIA PERCORRIDA EM \n     CRUZEIRO DA AERONAVE É \n        R = %E \n', R ) % [ dia ]


% Loiter-Fuel Fraction {RAYMER EQ.3.6} ==================|
Wloiter = exp( ( -E * SFC.FL250 ) / ( LDMAX ) ); % W4/W3 |
% =======================================================|

WxW0 = Wtakeoff * Wclimb * Wcruise * Wloiter *  Wlanding * Wdescent;

W.f0 = 1.05 * ( 1 - WxW0 ); % {RAYMER P.71}

% Takeof-Weight Guess
W0G = linspace( 2000, 20000, 100000 );

% Empty-Weight {RAYMER P.59}
W.e0 = 1.4 * W0G .^ ( -.1 ); 

% Takeof-Weight
W0  = ( Wcrew + Wpayload ) ./ ( 1 - W.f0 - W.e0 );

% DIFERENCE
WDif = W0-W0G;
[~,idx] = min(abs(W0-W0G));
WDif = WDif(1:idx);


%% GRÁFICOS

figure
plot(W0G(1:length(WDif)),W0(1:length(WDif)), 'LineWidth', 2)
hold on
plot(W0G(1:length(WDif)),WDif, '--','LineWidth', 1)
xline(W0G(length(WDif)),'-.','MAX WEIGHT')
xlabel('Takeoff-Weight Guess (We/W0) [-]')
ylabel('Takeoff-Weight (W0) [-]')
title('Takeoff-Weight Sizing')
legend('Takeoff-Weight',...
       'Difference',...
       'FontSize', 12,'FontName','Times New Roman','Location','best');
grid on; grid minor;
set(gcf, 'Color', 'w');
set(gca,'GridLineStyle', '-');
set(gcf,'paperPositionMode','auto')  

% DEFINE PESO
[~,idx] = min(abs(W0-W0G));
W0 = round(W0G(idx));

fprintf( '\n ============================ \n')
fprintf( '\n     O peso da aeronave é \n       W = %E \n', W0 ) % [ dia ]
fprintf('\n ======== FIM %s ============= \n',char( double( '👍😎' ) ) )
% =======================================================================
