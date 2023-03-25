%% Inicialização
% Referência: Gundmunssen
clear;clc;close all;
 
%% Parâmetros iniciais propostos para a aeronave

% Constantes
g = 9.81 ;       % Aceleração da gravidade [m/s²]

% Altitudes de Operação
h = 0;           % Altitude nível do mar [m]
hc = 45000.3048; % Teto operacional da aeronave [m] {ESTIMADO DADOS DA TABELA}

% Velocidades de Operação
Vc = 850 / 3.6; % Velocidade média de cruzeiro da aeronave [m/s] {ESTIMADO DADOS DE CONCORRENTES}
Vv = 20.80;   % Taxa de subida [m/s] {ESTIMADO DADOS DA TABELA}

%% Parâmetros iniciais propostos para a aeornave
% Geometria da Aeronave
AR   = 9;                                                           % Razão de Aspecto {ESTIMADO DADOS DE CONCORRENTES}
L_le = 28;                                                          % Angulo de Enflechamento {GRÁFICO RAYMER P.69}
e.sweep0   = 1.78 * ( 1 - .045 * (AR ^ .68) ) - .64;                % Eficiência de Evergadura asa reta
e.sweepvar = 4.61 * (1 - .045 * (AR ^ .68)) * (cosd(30)^.15) - 3.1; % Eficiência de Evergadura
e.interp   = interp1([0 30],[.7831 .5071],[28],'linear');           % Interpolação para melhor obtenção de eficiencia de Evergadura

%% Variação de W/S e T/W para diagrama de Gudmundsen

% Parâmetro de carga alar
GU.WS = linspace( 500, 5000, 10000 ); % Vetor Carga Alar [N/m²] {ESCALA BASEADA NO GUDMUNSON P.70}

% Distancias de Pista de Operação
Xd = 950; % Distância mínima de decolagem [m] {ESTIMADO DADOS DA TABELA}

% Parâmetros de Pista
mu = .04; % Coeficiente de atrito da pista {GUM TABLE 17-3}

% Atmosfera Padrão
[T, a, P, rho] = atmosisa([h hc]);

% Derivadas da Aeronave
      Clmax = 2.7;           % Coeficiente de sustentação máximo (Aula 06 - pag. 42)
 
    % Arrasto
      Cfe    = .003;                    % Coeficiente de Fricçào de Superfície Equivalente - (04 Estimativa de Polar de Arrasto)
      Cd_0   = .02;                      % Coeficiente de arrasto parasita ------------------ {DADO TABELADO MOHAMMAD P.70}
      k      = 1 / (e.interp * pi * AR); % Coeficiente de arrasto induzido ------------------ (04 Estimativa de Polar de Arrasto)
      Cd_min = 2 * Cd_0;                 % Coeficiente de arrasto mínimo -------------------- (04 Estimativa de Polar de Arrasto)

    % Decolagem
    % Cl = (2 .* GU.W)/(rho(2) * Sw * Vc); % Cl de decolagem 
      Cl = Clmax/1.44;                     % Cl de decolagem 
      Cd = Cd_0 + k .* Cl^2;               % Cd de decolagem                                            

% Requisito de Desempenho
Ps = 9.144; % Energia Excedente {POR NORMA 500 ft/min} [ m/s ] 

% Velocidades e Pressão dinâmica de Operação
qc = 0.5 * Vc ^2 * rho(2);                          % Pressão dinâmica (FLIGHT LEVEL)
Vs = sqrt( ( (2 / ( qc * Clmax ) ) .* GU.WS ) ); % Velocidade de Stall da aeronave [m/s]
Vdec = 1.2 .* Vs;                              % Velocidade máxima de decolagem [m/s]
qdec = .5 * Vdec .^ 2 * rho(1);                     % Pressão dinâmica (TAKEOFF)

%% Equações de Gudmundsen

% T/W for a Level Constant-velocity Turn - Pág.68 (Eq.3-1)

GU.TW_cvt = qc * ( Cd_min ./ ( GU.WS ) + ...
    k *( ( 1 / cosd(45) ) / qc )^2 .* ( GU.WS ) );

% T/W for a Desired Specific Energy Level - Pág. 68 (Eq. 3-2)

GU.TW_sel = qc * ( Cd_min ./ ( GU.WS ) + k *( ...
    ( 1/cosd(45) ) / qc )^2 .* ( GU.WS ) ) + Ps/Vc;

% T/W for a Desired Rate of Climb - pág.69 (Eq.3-3)

GU.TW_rc = Vv/Vc + ( ( qc * Cd_min ) ./ GU.WS ) + ( ( k .* GU.WS ) / qc ); 

% T/W for a Desired T-O Distance - pág.69 (Eq.3-4)

GU.TW_to = ( Vdec.^2 ./( 2 * g * Xd ) ) + ( qdec .* Cd ./ ( GU.WS ) )+ ...
           ( mu * ( 1 - qdec .* Cl ./ ( GU.WS ) ) );
          
% T/W for a Desired Cruise Airspeed - pág.59 (Eq.3-5)

GU.TW_cr = qc * Cd_min * ( 1 ./ ( GU.WS ) ) + ...
           k * ( 1 / qc ) .* ( GU.WS );
           
% T/W for a Service Ceiling (ROC 100 fpm or 0.508 m/s) - pág.69 (Eq.3-6)

GU.TW_ts = Vv ./ sqrt( ( ( 2 ./ rho(2) ) .* ( GU.WS ) .* sqrt( k / (3 * Cd_min) ) ) ) +...
           4 * sqrt( (k * Cd_min) / 3);

%% Plotando o Diagrama de Gudmundsen

figure
hold on
grid on
plot(GU.WS,GU.TW_cvt)
plot(GU.WS,GU.TW_sel)
plot(GU.WS,GU.TW_rc)
plot(GU.WS,GU.TW_to)
plot(GU.WS,GU.TW_cr)
plot(GU.WS,GU.TW_ts)
xlabel('Carga Alar - W/S [N/m²]')
ylabel('Empuxo Específico - T/W [-]')
title('DIAGRAMA DE GUNDMUNDSEN')
legend('T/W for a Level Constant-velocity Turn',...
       'T/W for a Desired Specific Energy Level',...
       'T/W for a Desired Rate of Climb',...
       'T/W for a Desired T-O Distance',...
       'T/W for a Desired Cruise Airspeed',...
       'T/W for a Service Ceiling (ROC 100 fpm or 0.508 m/s)', ...
       'Location','NorthWest','FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w');
set(gca,'GridLineStyle', '-');
set(gcf,'paperPositionMode','auto')  
ylim([0 1])
xline(2903,'-.','REQUISITO')
yline(0.2511,'-.','REQUISITO')