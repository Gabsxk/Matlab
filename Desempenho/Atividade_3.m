%% Inicialização 🤔🛠️
clc; close all; clear all; format compact;
% =======================================================================

fprintf( '\n ============================================= \n')
fprintf( '\n Atividade 1 de entrega \n')
fprintf( '\n Item b. FL330 \n')
fprintf( '\n ============================================= \n')
%% Atividade de sala 📚 @1
% Considerando a aeronave Embraer 170, propulsionada por dois motores
%  turbofan General Electric CF34-8E5, em condições de voo de cruzeiro
%  típico, com as especificações apresentadas abaixo, determinar as
%  velocidades polares máxima e mínima, a velocidade de arrasto total
%  mínimo e o excesso de empuxo máximo.

% Dados
% - Polar de arrasto da aeronave: CD = 0.0184 + 0.0346 * CL^2
% - Condição ambiente de cruzeiro: FL330 (33.000 pés), 0.409 kg/m^3
% - Área de asa: 72,7 m^2
% - Peso máximo de decolagem (MTOW): 353062 N
% - Peso médio em voo de cruzeiro (W): 307274 N
% - Empuxo máximo dos motores (ISA-SL @ 100 %TR): 128600 N
% - Empuxo dos motores em cruzeiro (FL330 @ 85 %TR): 22448 N

% Expressão do arrasto total em termos de CL
% - k * W * CL^2 - T * CL + CD0 * W = 0

% Equacionamento para condição de voo de cruzeiro
syms CL;
k = 0.0346;   % [ - ] Coeficiente de arrasto induzido da aeronave
W = 307274;   % [ N ] Peso médio em voo de cruzeiro 
T = 22448;    % [ N ] Empuxo dos motores em cruzeiro
CD0 = 0.0184; % [ - ] Coeficiente de arrasto parasita da aeronave

% Obtenção do coeficente de sustentação para condição de cruzeiro
eqn = k * W * CL^2 - T * CL + CD0 * W == 0; % Equação do arrasto total
CL = solve(eqn, CL);                        % Solução da equação
CL = double(CL);                            % Conversão para double
fprintf('\n Coeficiente de Sustentação da Aeronava para velocidade máxima: CL = %E [ - ] \n', CL(1)); % Saída do resultado
fprintf('\n Coeficiente de Sustentação da Aeronava para velocidade mínima: CL = %E [ - ] \n', CL(2)); % Saída do resultado

% Obtenção das velocidades máxima e mínima para condição de cruzeiro
rho = 0.409;                            % [kg/m^3] Densidade do ar
S = 72.7;                               % [ m^2 ] Área da asa
Vmax = sqrt(2 * W / (rho * S * CL(1))); % [ m/s ]Velocidade máxima
Vmin = sqrt(2 * W / (rho * S * CL(2))); % [ m/s ] Velocidade mínima
fprintf('\n Velocidade máxima da aeronave: Vmax = %E [ m/s ] \n', Vmax);   % [ m/s ] Saída do resultado de velocidade
fprintf('\n Velocidade mínima da aeronave: Vmin = %E [ m/s ] \n', Vmin);   % [ m/s ] Saída do resultado de velocidade

% Obtenção do minimo coeficente de sustentação de mínimo arrasto e arrasto mínimo
CLD0min = sqrt( ( CD0 ) / ( k ) );      % Coeficiente de sustentação para arrasto mínimo
CDmin = CD0 + k * CLD0min^2;            % Coeficiente de arrasto mínimo
fprintf('\n Coeficiente de sustentação de arrasto total mínimo: CLD0min = %E [ - ]\n', CLD0min); % Saída do resultado
fprintf('\n Coeficiente de arrasto total mínimo: CDmin = %E [ - ]\n', CDmin);                      % Saída do resultado

% Obtenção da velocidade para arrasto mínimo e força de arrasto mínimo
VminCD = sqrt( 2 * W / (rho * S * CLD0min)); % [ m/s ] Velocidade para arrasto mínimo
D = 0.5 * rho * S * VminCD^2 * CDmin;        % [ N   ] Força de arrasto mínima
fprintf('\n Velocidade de arrasto total mínimo: VminCD = %E [ m/s ] \n', VminCD);    % [ m/s ] Saída do resultado
fprintf('\n Arrasto total mínimo: D = %E [ N ] \n', D);                % [ N   ] Saída do resultado

% Obtenção do excesso de empuxo máximo
EEM = T - D;                      % [ N ] Excesso de empuxo máximo
fprintf('\n Excesso de empuxo máximo: DT = %E [ N ] \n', EEM); % [ N ] Saída do resultado

% Obtenção da potência de arrasto minima
CLPmin = sqrt( ( 3*CD0 ) / ( k ) );                 % [ -   ] Coeficiente de sustentação para potência de arrasto mínimo
CDPmin = CD0 + k * CLPmin^2;                        % [ -   ] Coeficiente de arrasto para potência de arrasto mínimo
VPmin = sqrt( 2 * W / ( rho * S * CLPmin ));        % [ m/s ] Velocidade para potência de arrasto mínimo
DPmin = 0.5 * rho * S * VPmin^2 * CDPmin;           % [ N   ] Força de arrasto para potência de arrasto mínima 
EEMP  = T * VPmin - 0.5 * rho * S * VPmin^3 * CD0;  % [ N   ] Excesso de empuxo máximo para potência de arrasto mínimo
fprintf('\n Velocidade para potência de arrasto mínima: VPmin = %E [ m/s ] \n', VPmin);      % [ N   ] Saída do resultado
fprintf('\n Força de arrasto para potência de arrasto mínima : DPmin = %E [ N ] \n', DPmin);     % [ N   ] Saída do resultado
fprintf('\n Excesso de empuxo máximo para potência de arrasto mínimo: DTP = %E [ N ] \n', EEMP); % [ N   ] Saída do resultado 

% Comportamento Vetorial
V    = linspace(.0001, 1000, 10000);        % [ m/s ] Vetor de velocidades 
CL   = 2 * W ./ ( rho * S * V.^2 );         % [ -   ]Coeficiente de sustentação
CD   = CD0 + k * CL .^2 ;                   % [ -   ]Coeficiente de arrasto
D   = (0.5 * rho * S * V.^2 .* CD)  * 1E-3; % [ kN  ] Força de arrasto total
D2  = (CD0 * W ./ CL + k * W .* CL) * 1E-3; % [ kN  ] Força de arrasto total
P    = T  .* V * 1E-3;                      % [ kW  ] Potência de empuxo
PD   = D .* V;                              % [ kW  ] Potência de arrasto
PD2  = D2 .* V;                             % [ kW  ] Potência de arrasto
EEM  = T * 1E-3 - D;                        % [ kN  ] Excesso de empuxo máximo
EEM2 = T * 1E-3 - D2;                       % [ kN  ] Excesso de empuxo máximo
EEMP = (T .* V - 0.5 * rho * S * V .^3 .*  CD0) * 1E-3;  % [ kN   ] Excesso de empuxo máximo para potência de arrasto mínimo

% Plotagem
figure; 
subplot(2,2,1);grid minor; hold on;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
plot(V * 3.6, D2, 'k--', 'LineWidth', 2 );
text(V(600) * 3.6, D(600),'\leftarrow Total Drag (D)')
yline(T * 1E-3, 'b-.', 'Thrust Avalible', 'LineWidth', 1 );
ylim( [ 0 1E2 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Thrust or Total Drag [kN]')
title('Polar de Arrasto ')
legend( 'Total Drag by V',...
        'Total Drag by CL',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,2);grid minor; hold on;
plot(V * 3.6, P, 'b--','LineWidth',  1 );
text(V(500) * 3.6, P(500),'\leftarrow Propulsive Power (TV)')
plot(V * 3.6, PD, 'r-','LineWidth', .5 );
plot(V * 3.6, PD2, 'k--','LineWidth', 2 );
text(V(300) * 3.6, PD2(300),'\leftarrow Drag Power (DV)')
ylim( [ 0 1E4 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Propulsive or Drag Power [kN]')
title(' Curva de Potência de Polar de Arrasto ')
legend( 'Propulsive Power',...
        'Drag Power by V',...
        'Drag Power by CL',...
        'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,3);grid minor; hold on;
yyaxis left;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
plot(V * 3.6, D2, 'k--', 'LineWidth', 2 );
text(V(500) * 3.6, D(500),'\leftarrow Total Drag (D)')
yline(T * 1E-3, 'b-.', 'Thrust Avalible', 'LineWidth', 1 );
ylim( [ 0 3E2 ] );
ylabel('Thrust or Total Drag [kN]')
yyaxis right;
plot(V * 3.6, EEM, 'g-','LineWidth', .5 );
plot(V * 3.6, EEM2, 'k--', 'LineWidth', 2 );
text(V(1800) * 3.6, EEM(1800),'\leftarrow Excess Propulsive Thrust')
ylim( [ 0 .1E2 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Excess Propulsive Thrust [kN]')
title('Polar de Arrasto ')
legend( 'Total Drag by V',...
        'Total Drag by CL',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,4);grid minor; hold on;
yyaxis left;
plot(V * 3.6, P, 'b--','LineWidth',  1 );
text(V(200) * 3.6, P(200),'\leftarrow Propulsive Power (TV)')
plot(V * 3.6, PD, 'r-','LineWidth', .5 );
plot(V * 3.6, PD2, 'k--','LineWidth', 2 );
text(V(300) * 3.6, PD2(300),'\leftarrow Drag Power (DV)')
ylabel('Propulsive or Drag Power [kN]')
ylim( [ 0 1E4 ] );
yyaxis right;
plot(V * 3.6, EEMP, 'g--','LineWidth',  1 );
text(V(2500) * 3.6, EEMP(2500),'\leftarrow Excess Propulsive Power')
ylabel('Excess Propulsive Power [kN]')
ylim( [ 0 .5E4 ] );
xlabel('True Airspeed [Km/h]')
title(' Curva de Potência de Polar de Arrasto ')
legend( 'Propulsive Power',...
        'Drag Power by V',...
        'Drag Power by CL',...
        'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

% =======================================================================

fprintf( '\n ============================================= \n')
fprintf( '\n            Atividade 2 de sala \n')
fprintf( '\n                    Item b. \n')
fprintf( '\n ============================================= \n')

%% Atividade de sala 📚 @2
% Considerando a aeronave Embraer BEM-121 Xingu, propulsionada pelo
%  motor turbohélice PW&C PT6A-135, em condições de voo de cruzeiro
%  típico, determinar as velocidades polares mínima e máxima, a
%  velocidade de potência de arrasto mínima, a velocidade de excesso de
%  potência máximo e o excesso de potência.

% Dados
% - Polar de arrasto da aeronave: CD = 0.0276 + 0.0441 * CL^2
% - Condição ambiente de cruzeiro: FL200 (20.000 pés), rho = 0.652 kg/m^3
% - Área de asa: 27,5 m^2
% - Peso máximo de decolagem (MTOW): 55623 N
% - Peso médio em voo de cruzeiro (W): 48839 N
% - Empuxo máximo dos motores (ISA-SL @ 100 %TR): 1500 shp
% - Empuxo dos motores em cruzeiro (357172W @ 75 %TR): 357172 W

% Expressão do arrasto total em termos de V
% - .5 * rho * V^4 * S * CD0 - V * tep + 2 * k * W ^2 / ( rho * S )  == 0

% Equacionamento para condição de voo de cruzeiro
syms V;
S   = 27.5;         % [ m^2    ] Área de asa da aeronave
CD0 = .0276;        % [ -      ] Coeficiente de arrasto parasita da aeronave
tep = 1500 * 745.7; % [ W      ] Empuxo dos motores em cruzeiro
k   = .0441;        % [ -      ] Coeficiente de arrasto induzido da aeronave
W   = 48839;        % [ N      ] Peso médio em voo de cruzeiro 

% Condições atmosféricas FL100
[~, ~, ~, rho] = atmosisa( 10000 * .3048) 
tep = tep * rho / 1.225 ;

% Obtenção do coeficente de sustentação para condição de cruzeiro
eqn = .5 * rho * V^4 * S * CD0 - V * tep + 2 * k * W ^2 / ( rho * S ) == 0; % Equação do arrasto total
V = solve(eqn, V);                        % Solução da equação
V = double(V);                            % Conversão para double
fprintf('\n Velocidade máxima de voo da aeornave: V = %E [ m/s ] \n', max( V(V>=0) ) );           % Saída do resultado
fprintf('\n Velocidade mínima de voo da aeornave: V = %E [ m/s ] \n', min( V(V>=0) ) );           % Saída do resultado

% Obtenção do minimo coeficente de sustentação de mínimo arrasto e arrasto mínimo
CLD0min = sqrt( ( CD0 ) / ( k ) );      % Coeficiente de sustentação para arrasto mínimo
CDmin = CD0 + k * CLD0min^2;            % Coeficiente de arrasto mínimo
fprintf('\n Coeficiente de sustentação de arrasto total mínimo: CLD0min = %E [ - ]\n', CLD0min); % Saída do resultado
fprintf('\n Coeficiente de arrasto total mínimo: CDmin = %E [ - ]\n', CDmin);                      % Saída do resultado

% Obtenção da velocidade para arrasto mínimo e força de arrasto mínimo
VminCD = sqrt( 2 * W / (rho * S * CLD0min)); % [ m/s ] Velocidade para arrasto mínimo
D = 0.5 * rho * S * VminCD^2 * CDmin;        % [ N   ] Força de arrasto mínima
fprintf('\n Velocidade de arrasto total mínimo: VminCD = %E [ m/s ] \n', VminCD);    % [ m/s ] Saída do resultado
fprintf('\n Arrasto total mínimo: D = %E [ N ] \n', D);                % [ N   ] Saída do resultado

% Obtenção de propriedades para potência mínima de arrasto
CLPmin = sqrt( ( 3 * CD0 ) / ( k ) );          % [ -   ] Coeficiente de sustentação para potência de arrasto mínimo
CDPmin = CD0 + k * CLPmin^2;                   % [ -   ] Coeficiente de arrasto para potência de arrasto mínimo
VPmin  = sqrt( 2 * W / ( rho * S * CLPmin ) ); % [ m/s ] Velocidade para potência de arrasto mínimo
DPmin  = .5 * rho * VPmin^2 * S * CDPmin;      % [ N   ] Força de arrasto para potência de arrasto mínima 
PDmin  = VPmin * DPmin;                        % [ W   ] Potência de arrasto mínimo
EEM = tep - PDmin;                             % [ W   ] Excesso de empuxo máximo
fprintf('\n Coeficiente de sustentação para potência de arrasto mínimo: CLPmin = %E [ - ] \n', CLPmin); % [ m/s ] Saída do resultado 
fprintf('\n Velocidade para potência de arrasto mínima: VPmin = %E [ m/s ] \n', VPmin);                 % [ m/s   ] Saída do resultado
fprintf('\n Força de arrasto para potência de arrasto mínima : DPmin = %E [ N ] \n', DPmin);            % [ N   ] Saída do resultado
fprintf('\n Potência de arrasto mínimo: PDmin = %E [ N ] \n', PDmin);                                     % [ W   ] Saída do resultado 
fprintf('\n Excesso de empuxo máximo: Dtep = %E [ W ] \n', EEM);                                        % [ N   ] Saída do resultado

% Comportamento Vetorial
V    = linspace(.0001, 1000, 10000);           % [ m/s ] Vetor de velocidades 
CL   = 2 * W ./ ( rho * S * V.^2 );            % [ -   ] Coeficiente de sustentação
CD   = CD0 + k * CL .^2 ;                      % [ -   ] Coeficiente de arrasto
D    = (0.5 * rho * S .* V .^2 .* CD)  * 1E-3; % [ kN  ] Força de arrasto total
T    = tep * 1E-3 ./ V ;                       % [ kN  ] Empuxo total
EEM  = T - D ;                                 % [ kN  ] Excesso de empuxo máximo

% Plotagem
figure;grid minor; hold on;
yyaxis left;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
text(V(600) * 3.6, D(600),'\leftarrow Total Drag (D)')
plot(V * 3.6, T, 'b-', 'LineWidth', 1 );
text(V(600) * 3.6, P(600),'\leftarrow Thrust Avalible (T)')
ylim( [ 0 4E+1 ] )
ylabel('Thrust or Total Drag [kN]')
yyaxis right;
plot( V * 3.6, EEM, 'k-','LineWidth', 2 );
text(V(600) * 3.6, EEM(600),'\leftarrow Excess Propulsive Thrust')
ylim( [ 0 3E+1 ] )
xlim( [ 0 9E+2 ] )
xlabel('True Airspeed [Km/h]')
ylabel('Excess Propulsive Thrust [kN]')
title('Polar de Empuxo ')
legend( 'Total Drag by V',...
        'Thrust Avalibe',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

% =======================================================================
fprintf( '\n ============================================= \n')
fprintf( '\n            Atividade 1 de Entrega \n')
fprintf( '\n                Item a. FL0 \n')
fprintf( '\n ============================================= \n')

%% Atividade de sala 📚 @1
%  Utilizando os dados da aeronave a jato do Exemplo 2 (slide 41), 
%    plote as curvas polares de empuxo e potência para a condição 
%    SL.

% Dados
% - Polar de arrasto da aeronave: CD = 0.0184 + 0.0346 * CL^2
% - Condição ambiente de cruzeiro: FL330 (33.000 pés), 0.409 kg/m^3
% - Área de asa: 72,7 m^2
% - Peso máximo de decolagem (MTOW): 353062 N
% - Peso médio em voo de cruzeiro (W): 307274 N
% - Empuxo máximo dos motores (ISA-SL @ 100 %TR): 128600 N
% - Empuxo dos motores em cruzeiro (FL330 @ 85 %TR): 22448 N

% Expressão do arrasto total em termos de CL
% - k * W * CL^2 - T * CL + CD0 * W = 0

% Equacionamento para condição de voo de cruzeiro
syms CL;
k = 0.0346;   % [ - ] Coeficiente de arrasto induzido da aeronave
W = 307274;   % [ N ] Peso médio em voo de cruzeiro 
T = 128600;    % [ N ] Empuxo dos motores em cruzeiro para FL-SL
CD0 = 0.0184; % [ - ] Coeficiente de arrasto parasita da aeronave

% Obtenção das propriedades do ar para ISA-SL
[ ~, ~, ~, rho ] = atmosisa( 0 ); 

% Obtenção do coeficente de sustentação para condição de cruzeiro
eqn = k * W * CL^2 - T * CL + CD0 * W == 0; % Equação do arrasto total
CL = solve(eqn, CL);                        % Solução da equação
CL = double(CL);                            % Conversão para double
fprintf('\n Coeficiente de Sustentação da Aeronava para velocidade máxima: CL = %E [ - ] \n', CL(1)); % Saída do resultado
fprintf('\n Coeficiente de Sustentação da Aeronava para velocidade mínima: CL = %E [ - ] \n', CL(2)); % Saída do resultado

% Obtenção das velocidades máxima e mínima para condição de cruzeiro
S = 72.7;                               % [ m^2 ] Área da asa
Vmax = sqrt(2 * W / (rho * S * CL(1))); % [ m/s ]Velocidade máxima
Vmin = sqrt(2 * W / (rho * S * CL(2))); % [ m/s ] Velocidade mínima
fprintf('\n Velocidade máxima da aeronave: Vmax = %E [ m/s ] \n', Vmax);   % [ m/s ] Saída do resultado de velocidade
fprintf('\n Velocidade mínima da aeronave: Vmin = %E [ m/s ] \n', Vmin);   % [ m/s ] Saída do resultado de velocidade

% Obtenção do minimo coeficente de sustentação de mínimo arrasto e arrasto mínimo
CLD0min = sqrt( ( CD0 ) / ( k ) );      % Coeficiente de sustentação para arrasto mínimo
CDmin = CD0 + k * CLD0min^2;            % Coeficiente de arrasto mínimo
fprintf('\n Coeficiente de sustentação de arrasto total mínimo: CLD0min = %E \n', CLD0min); % Saída do resultado
fprintf('\n Coeficiente de arrasto total mínimo: CDmin = %E \n', CDmin);                      % Saída do resultado

% Obtenção da velocidade para arrasto mínimo e força de arrasto mínimo
VminCD = sqrt( 2 * W / (rho * S * CLD0min)); % [ m/s ] Velocidade para arrasto mínimo
D = 0.5 * rho * S * VminCD^2 * CDmin;        % [ N   ] Força de arrasto mínima
fprintf('\n Velocidade de arrasto total mínimo: VminCD = %E [ m/s ] \n', VminCD);    % [ m/s ] Saída do resultado
fprintf('\n Arrasto total mínimo: D = %E [ N ] \n', D);                % [ N   ] Saída do resultado

% Obtenção do excesso de empuxo máximo
EEM = T - D;                      % [ N ] Excesso de empuxo máximo
fprintf('\n Excesso de empuxo máximo: DT = %E N \n', EEM); % [ N ] Saída do resultado

% Obtenção da potência de arrasto minima
CLPmin = sqrt( ( 3*CD0 ) / ( k ) );                 % [ -   ] Coeficiente de sustentação para potência de arrasto mínimo
CDPmin = CD0 + k * CLPmin^2;                        % [ -   ] Coeficiente de arrasto para potência de arrasto mínimo
VPmin = sqrt( 2 * W / ( rho * S * CLPmin ));        % [ m/s ] Velocidade para potência de arrasto mínimo
DPmin = 0.5 * rho * S * VPmin^2 * CDPmin;           % [ N   ] Força de arrasto para potência de arrasto mínima 
EEMP  = T * VPmin - 0.5 * rho * S * VPmin^3 * CD0;  % [ N   ] Excesso de empuxo máximo para potência de arrasto mínimo
fprintf('\n Velocidade para potência de arrasto mínima: VPmin = %E [ m/s ] \n', VPmin);      % [ N   ] Saída do resultado
fprintf('\n Força de arrasto para potência de arrasto mínima : DPmin = %E [ N ] \n', DPmin);     % [ N   ] Saída do resultado
fprintf('\n Excesso de empuxo máximo para potência de arrasto mínimo: DTP = %E [ N ] \n', EEMP); % [ N   ] Saída do resultado 

% Comportamento Vetorial
V    = linspace(.0001, 1000, 10000);        % [ m/s ] Vetor de velocidades 
CL   = 2 * W ./ ( rho * S * V.^2 );         % [ -   ]Coeficiente de sustentação
CD   = CD0 + k * CL .^2 ;                   % [ -   ]Coeficiente de arrasto
D   = (0.5 * rho * S * V.^2 .* CD)  * 1E-3; % [ kN  ] Força de arrasto total
D2  = (CD0 * W ./ CL + k * W .* CL) * 1E-3; % [ kN  ] Força de arrasto total
P    = T  .* V * 1E-3;                      % [ kW  ] Potência de empuxo
PD   = D .* V;                              % [ kW  ] Potência de arrasto
PD2  = D2 .* V;                             % [ kW  ] Potência de arrasto
EEM  = T * 1E-3 - D;                        % [ kN  ] Excesso de empuxo máximo
EEM2 = T * 1E-3 - D2;                       % [ kN  ] Excesso de empuxo máximo
EEMP = (T .* V - 0.5 * rho * S * V .^3 .*  CD0) * 1E-3;  % [ kN   ] Excesso de empuxo máximo para potência de arrasto mínimo

% Plotagem
figure; 
subplot(2,2,1);grid minor; hold on;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
plot(V * 3.6, D2, 'k--', 'LineWidth', 2 );
text(V(600) * 3.6, D(600),'\leftarrow Total Drag (D)')
yline(T * 1E-3, 'b-.', 'Thrust Avalible', 'LineWidth', 1 );
ylim( [ 0 3E2 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Thrust or Total Drag [kN]')
title('Polar de Arrasto ')
legend( 'Total Drag by V',...
        'Total Drag by CL',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,2);grid minor; hold on;
plot(V * 3.6, P, 'b--','LineWidth',  1 );
text(V(500) * 3.6, P(500),'\leftarrow Propulsive Power (TV)')
plot(V * 3.6, PD, 'r-','LineWidth', .5 );
plot(V * 3.6, PD2, 'k--','LineWidth', 2 );
text(V(300) * 3.6, PD2(300),'\leftarrow Drag Power (DV)')
ylim( [ 0 1E4 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Propulsive or Drag Power [kN]')
title(' Curva de Potência de Polar de Arrasto ')
legend( 'Propulsive Power',...
        'Drag Power by V',...
        'Drag Power by CL',...
        'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,3);grid minor; hold on;
yyaxis left;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
plot(V * 3.6, D2, 'k--', 'LineWidth', 2 );
text(V(500) * 3.6, D(500),'\leftarrow Total Drag (D)')
yline(T * 1E-3, 'b-.', 'Thrust Avalible', 'LineWidth', 1 );
ylim( [ 0 3E2 ] );
ylabel('Thrust or Total Drag [kN]')
yyaxis right;
plot(V * 3.6, EEM, 'g-','LineWidth', .5 );
plot(V * 3.6, EEM2, 'k--', 'LineWidth', 2 );
text(V(1800) * 3.6, EEM(1800),'\leftarrow Excess Propulsive Thrust')
ylim( [ 0 1.5E2 ] );
xlabel('True Airspeed [Km/h]')
ylabel('Excess Propulsive Thrust [kN]')
title('Polar de Arrasto ')
legend( 'Total Drag by V',...
        'Total Drag by CL',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

subplot(2,2,4);grid minor; hold on;
yyaxis left;
plot(V * 3.6, P, 'b--','LineWidth',  1 );
text(V(200) * 3.6, P(200),'\leftarrow Propulsive Power (TV)')
plot(V * 3.6, PD, 'r-','LineWidth', .5 );
plot(V * 3.6, PD2, 'k--','LineWidth', 2 );
text(V(300) * 3.6, PD2(300),'\leftarrow Drag Power (DV)')
ylabel('Propulsive or Drag Power [kN]')
ylim( [ 0 1E4 ] );
yyaxis right;
plot(V * 3.6, EEMP, 'g--','LineWidth',  1 );
text(V(2500) * 3.6, EEMP(2500),'\leftarrow Excess Propulsive Power')
ylabel('Excess Propulsive Power [kN]')
ylim( [ 0 2.5E4 ] );
xlabel('True Airspeed [Km/h]')
title(' Curva de Potência de Polar de Arrasto ')
legend( 'Propulsive Power',...
        'Drag Power by V',...
        'Drag Power by CL',...
        'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

% =======================================================================

fprintf( '\n ============================================= \n')
fprintf( '\n            Atividade 2 de sala \n')
fprintf( '\n                Item a. FL0 \n')
fprintf( '\n ============================================= \n')

%% Atividade de sala 📚 @2
% Considerando a aeronave Embraer BEM-121 Xingu, propulsionada pelo
%  motor turbohélice PW&C PT6A-135, em condições de voo de cruzeiro
%  típico, determinar as velocidades polares mínima e máxima, a
%  velocidade de potência de arrasto mínima, a velocidade de excesso de
%  potência máximo e o excesso de potência.

% Dados
% - Polar de arrasto da aeronave: CD = 0.0276 + 0.0441 * CL^2
% - Condição ambiente de cruzeiro: FL200 (20.000 pés), rho = 0.652 kg/m^3
% - Área de asa: 27,5 m^2
% - Peso máximo de decolagem (MTOW): 55623 N
% - Peso médio em voo de cruzeiro (W): 48839 N
% - Empuxo máximo dos motores (ISA-SL @ 100 %TR): 1500 shp
% - Empuxo dos motores em cruzeiro (357172W @ 75 %TR): 357172 W

% Expressão do arrasto total em termos de V
% - .5 * rho * V^4 * S * CD0 - V * tep + 2 * k * W ^2 / ( rho * S )  == 0

% Equacionamento para condição de voo de cruzeiro
syms V;
S   = 27.5;         % [ m^2    ] Área de asa da aeronave
CD0 = .0276;        % [ -      ] Coeficiente de arrasto parasita da aeronave
tep = 1500 * 745.7; % [ W      ] Empuxo máximo dos motores
k   = .0441;        % [ -      ] Coeficiente de arrasto induzido da aeronave
W   = 48839;        % [ N      ] Peso médio em voo de cruzeiro 

% Obtenção das propriedades do ar para ISA-SL
[ ~, ~, ~, rho ] = atmosisa( 0 ); 

% Obtenção do coeficente de sustentação para condição de cruzeiro
eqn = .5 * rho * V^4 * S * CD0 - V * tep + 2 * k * W ^2 / ( rho * S ) == 0; % Equação do arrasto total
V = solve(eqn, V);                        % Solução da equação
V = double(V);                            % Conversão para double
fprintf('\n Velocidade máxima de voo da aeornave: V = %E m/s \n', max( V(V>=0) ) );           % Saída do resultado
fprintf('\n Velocidade mínima de voo da aeornave: V = %E m/s \n', min( V(V>=0) ) );           % Saída do resultado

% Obtenção do minimo coeficente de sustentação de mínimo arrasto e arrasto mínimo
CLD0min = sqrt( ( CD0 ) / ( k ) );      % Coeficiente de sustentação para arrasto mínimo
CDmin = CD0 + k * CLD0min^2;            % Coeficiente de arrasto mínimo
fprintf('\n Coeficiente de sustentação de arrasto total mínimo: CLD0min = %E \n', CLD0min); % Saída do resultado
fprintf('\n Coeficiente de arrasto total mínimo: CDmin = %E \n', CDmin);                      % Saída do resultado

% Obtenção da velocidade para arrasto mínimo e força de arrasto mínimo
VminCD = sqrt( 2 * W / (rho * S * CLD0min)); % [ m/s ] Velocidade para arrasto mínimo
D = 0.5 * rho * S * VminCD^2 * CDmin;        % [ N   ] Força de arrasto mínima
fprintf('\n Velocidade de arrasto total mínimo: VminCD = %E [ m/s ] \n', VminCD);    % [ m/s ] Saída do resultado
fprintf('\n Arrasto total mínimo: D = %E [ N ] \n', D);                % [ N   ] Saída do resultado

% Obtenção de propriedades para potência mínima de arrasto
CLPmin = sqrt( ( 3 * CD0 ) / ( k ) );          % [ -   ] Coeficiente de sustentação para potência de arrasto mínimo
CDPmin = CD0 + k * CLPmin^2;                   % [ -   ] Coeficiente de arrasto para potência de arrasto mínimo
VPmin  = sqrt( 2 * W / ( rho * S * CLPmin ) ); % [ m/s ] Velocidade para potência de arrasto mínimo
DPmin  = .5 * rho * VPmin^2 * S * CDPmin;      % [ N   ] Força de arrasto para potência de arrasto mínima 
PDmin  = VPmin * DPmin;                        % [ W   ] Potência de arrasto mínimo
EEM = tep - PDmin;                             % [ W   ] Excesso de empuxo máximo
fprintf('\n Coeficiente de sustentação para potência de arrasto mínimo: CLPmin = %E [ - ] \n', CLPmin); % [ m/s ] Saída do resultado 
fprintf('\n Velocidade para potência de arrasto mínima: VPmin = %E [ m/s ] \n', VPmin);                 % [ m/s   ] Saída do resultado
fprintf('\n Força de arrasto para potência de arrasto mínima : DPmin = %E [ N ] \n', DPmin);            % [ N   ] Saída do resultado
fprintf('\n Potência de arrasto mínimo: PDmin = %E [ N ] \n', PDmin);                                     % [ W   ] Saída do resultado 
fprintf('\n Excesso de empuxo máximo: Dtep = %E [ W ] \n', EEM);                                        % [ N   ] Saída do resultado

% Comportamento Vetorial
V    = linspace(.0001, 1000, 1000000);           % [ m/s ] Vetor de velocidades 
CL   = 2 * W ./ ( rho * S * V.^2 );            % [ -   ] Coeficiente de sustentação
CD   = CD0 + k * CL .^2 ;                      % [ -   ] Coeficiente de arrasto
D    = (0.5 * rho * S .* V .^2 .* CD)  * 1E-3; % [ kN  ] Força de arrasto total
T    = tep * 1E-3 ./ V ;                       % [ kN  ] Empuxo total
EEM  = T - D ;                                 % [ kN  ] Excesso de empuxo máximo

% Plotagem
figure;grid minor; hold on;
yyaxis left;
plot(V * 3.6, D, 'r-','LineWidth', .5 );
text(V(600) * 3.6, D(600),'\leftarrow Total Drag (D)')
plot(V * 3.6, T, 'b-', 'LineWidth', 1 );
text(V(600) * 3.6, T(600),'\leftarrow Thrust Avalible (T)')
ylim( [ 0 3E+2 ] )
ylabel('Thrust or Total Drag [kN]')
yyaxis right;
plot( V * 3.6, EEM, 'k-','LineWidth', 2 );
text(V(600) * 3.6, EEM(600),'\leftarrow Excess Propulsive Thrust')
ylim( [ 0 6E+1 ] )
xlim( [ 0 5E2 ] )
xlabel('True Airspeed [Km/h]')
ylabel('Excess Propulsive Thrust [kN]')
title('Polar de Empuxo ')
legend( 'Total Drag by V',...
        'Thrust Avalibe',...
         'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')
fprintf('\n ====================== FIM %s \n',char( double( '👍😎' ) ) )