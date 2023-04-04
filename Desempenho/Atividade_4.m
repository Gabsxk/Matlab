%% Inicialização 🤔🛠️
clc; close all; clear all; format compact;
% =======================================================================

%% Atividade de sala 📚 @1 e @2
% 1. Determine que a máxima eficiência aerodinâmica
% 2. Determine que Emax acontece para a condição de mínima força de arrasto


% E = L/D = CL/CD = CL / ( CD_0 + K * CL^2 )
% 1/E = ( CD_0 / CL ) + ( K * CL )
% d(1/E)/dCL = (- CD_0 / CL^2) + ( K )
% (- CD_0 / CL^2) + ( K ) = 0 -> CD_0 = K * CL^2 -> CD_0 = CD_i
% CD_0 = K * CL^2 -> CL_md = sqrt( CD_0 / K )
% CD_min = CD_0 + K * CL_md^2 -> CD_min = 2 * CD_0

% E = L / D = CL / CD -> Emax = CL_md / CD_min -> Emax = sqrt( CD_0 / K ) / ( 2 * CD_0 )
% Emax = sqrt( 1 / K * CD_0^2 ) / ( 2  ) 
% Emax = 1 / ( 2 * sqrt( K * CD_0^2 ) )

% =======================================================================

%% Atividade de sala 📚 @2
% 1. Determine a velocidade e o CL para o máximo alcance específico de 
% uma aeronave em função da velocidade de força de arrasto mínima (V md).

% ( CL ^( 1 / 2 ) / CD  ) = ( CL ^( 1 / 2 ) / ( CD_0 + K * CL ^2 ) )
% ( CD_0 + K * CL ^2 ) / ( CL ^( 1 / 2 ) ) = ( CD_0  / ( CL ^( 1 / 2 ) ) + ( K * CL ^( 3 / 2 ) )
% d( ( CD_0  / ( CL ^( 1 / 2 ) ) + ( K * CL ^( 3 / 2 ) ) ) / d( CL ) = 0
% (- .5 * CD_0  / ( CL ^( 3 / 2 ) ) + ( ( 3 / 2 ) * K * CL ^( 1 / 2 ) ) = 0
% (  CD_0  ) = ( ( 3 ) * K * CL ^( 2 ) )
% (  CD_0 / ( 3 * K ) ) =  CL ^( 2 ) 
% CLmaxsar = sqrt( CD_0 / ( 3 * K ) )

% u = Vmind / Vmaxsar -> u = sqrt(CLmind / CLmaxsar) = sqrt( sqrt( CD_0 / K ) / sqrt( CD_0 / ( 3 * K ) ) )
% u =  sqrt( sqrt( CD_0 / K  / ( CD_0 / ( 3 * K ) ) ) -> u = ( 3 ) ^( 1 / 4 )

% =======================================================================

%% Atividade de sala 📚 @3
% l. Uma aeronave de transporte de passageiros a jato possui as
%     seguintes características:

%                    S = 341,5m2 CD = 0,016 + 0,065 * CL2; TSFC = 0.8 1/h

%    Para um voo a 30.000 ft, essa aeronave possui massa de 100.000 kg,
%     sendo 30.000 kg de massa de combustível.

%     a) Se a aeronave voa a uma velocidade de 325,8 knot, qual o
%         alcance e a autonomia que ela consegue atingir sem
%         reabastecimento?

%     b) Qual seria a velocidade que a aeronave deveria voar para obter o
%         máximo alcance e a máxima autonomia? E quais seriam os
%         valores do máximo alcance e da máxima autonomia?

% OBS: Considere velocidade e CL constantes.

% Dados 
g = 9.81;                  % [ m/s^2 ] Aceleração da gravidade
S = 341.5;                 % [ m^2 ] Área da asa
m = 100000;                % [ kg ] Massa da aeronave
m_F = 30000;               % [ kg ] Massa de combustível
h = 30000 * 3048;          % [ m ] Altitude
V = 325.8 * 0.514444;      % [ m/s ] Velocidade da aeornave
CD0 = 0.016;              % [ - ] Coeficiente de Arrasto
k = 0.065;                 % [ - ] Coeficiente de Arrasto
TSFC = .8 / 3600;          % [ 1/s ] Consumo específico de combustível 

% a) Alcance e autonomia

% Obtenção de dados de atmosfera
[~, ~, ~, rho] = atmosisa( h );

% Obtenção da velocidade e coeficiente de sustentação de mínimo arrasto
CLmd = sqrt( CD0 / k );
Vmd = sqrt( 2 * m * g  / ( CLmd * rho * S ) );

% Obtenção da Máxima Eficiencia Aerodinâmica
Emax = 1 / ( 2 * sqrt( k * CD0 ) );

% Obtenção da Razão de Velocidade
u = V / Vmd;

% Obtenção da Razão de peso de cruzeiro
ome = ( m * g ) / ( ( m - m_F ) * g );

% Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
R1.V = ( Vmd / TSFC * Emax) * ( 2 * u ^3 / ( u ^4 + 1 ) ) * log( ome );

% Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
E1.V = ( 1 / TSFC * Emax) * ( 2 * u ^2 / ( u ^4 + 1 ) ) * log( ome );

% Velocidade para máximo alcanse e Máximo Alcanse
Vmax(1) = Vmd * ( 3 ) ^ ( 1 / 4 );
u1R.max  = Vmax(1) / Vmd;
R1.Vmax = ( Vmd / TSFC * Emax) * ( 2 * u1R.max ^3 / ( u1R.max ^4 + 1 ) ) * log( ome );

% Velocidade para máxima autonomia e máxima Autonomia
Vmax(2) = Vmd * 1;
u1E.max  = Vmax(2) / Vmd;
E1.Vmax = ( 1 / TSFC * Emax) * ( 2 * u1E.max ^3 / ( u1E.max ^4 + 1 ) ) * log( ome );

% Ipressão de respostas 
fprintf('\n Atividade de Sala 3 \n');
fprintf('          Item a.      \n');
fprintf('\n Alcanse para velocidade de %E [ m/s ]: R = %E [ m ] \n', V, R1.V);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia para velocidade de %E [ m/s ]: E = %E [ s ] \n', V, E1.V);    % [ m/s ] Saída do resultado

fprintf('\n          Item b.      \n');
fprintf('\n Alcanse máximo para velocidade de máximo alcanse %E [ m/s ]: R = %E [ m ] \n', Vmax(1), R1.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia máxima para velocidade de máxima autonomia %E [ m/s ]: E = %E [ s ] \n', Vmax(2), E1.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n ============================================================== \n');    % [ m/s ] Saída do resultado

% ==========================================================================

%% Atividade de sala 📚 @4
% l. Uma aeronave de transporte de passageiros a jato possui as
%     seguintes características:

%                    S = 341,5m2 CD = 0,016 + 0,065CL2; TSFC = 0.8 1/h

%    Para um voo a 30.000 ft, essa aeronave possui massa de 100.000 kg,
%     sendo 30.000 kg de massa de combustível.

%     a) Se a aeronave voa a uma velocidade de 325,8 knot, qual o
%         alcance e a autonomia que ela consegue atingir sem
%         reabastecimento?

%     b) Qual seria a velocidade que a aeronave deveria voar para obter o
%         máximo alcance e a máxima autonomia? E quais seriam os
%         valores do máximo alcance e da máxima autonomia?

% OBS: Considere velocidade e CL constantes.

% Dados 
g = 9.81;                  % [ m/s^2 ] Aceleração da gravidade
S = 341.5;                 % [ m^2 ] Área da asa
m = 100000;                % [ kg ] Massa da aeronave
m_F = 30000;               % [ kg ] Massa de combustível
h = 30000 * 3048;          % [ m ] Altitude
V = 325.8 * 0.514444;      % [ m/s ] Velocidade da aeornave
CD0 = 0.016;              % [ - ] Coeficiente de Arrasto
k = 0.065;                 % [ - ] Coeficiente de Arrasto
TSFC = .8 / 3600;          % [ 1/s ] Consumo específico de combustível 

% a) Alcance e autonomia

% Obtenção de dados de atmosfera
[~, ~, ~, rho] = atmosisa( h );

% Obtenção da velocidade e coeficiente de sustentação de mínimo arrasto
CLmd = sqrt( CD0 / k );
Vmd = sqrt( 2 * m * g  / ( CLmd * rho * S ) );

% Obtenção da Máxima Eficiencia Aerodinâmica
Emax = 1 / ( 2 * sqrt( k * CD0 ) );

% Obtenção da Razão de Velocidade
u = V / Vmd;

% Obtenção da Razão de peso de cruzeiro
ome = ( m * g ) / ( ( m - m_F ) * g );

% Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
R2.V = ( Vmd / TSFC * Emax) * ( 2 * u ^3 / ( u ^4 + 1 ) ) * 2 * ( 1 - sqrt( ome ) ^( -1 ) );

% Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
E2.V = ( 1 / TSFC * Emax) * ( 2 * u ^2 / ( u ^4 + 1 ) ) * log( ome );

% Velocidade para máximo alcanse e Máximo Alcanse
Vmax(1) = Vmd * ( 3 ) ^ ( 1 / 4 );
u2R.max = Vmax(1) / Vmd;
R2.Vmax = ( Vmd / TSFC * Emax) * ( 2 * u2R.max ^3 / ( u2R.max ^4 + 1 ) ) * 2 * ( 1 - sqrt( ome ) ^( -1 ) );

% Velocidade para máxima autonomia e máxima Autonomia
Vmax(2) = Vmd * 1;
u2E.max  = Vmax(2) / Vmd;
E2.Vmax = ( 1 / TSFC * Emax) * ( 2 * u2E.max ^2 / ( u2E.max ^4 + 1 ) ) * log( ome );

% Ipressão de respostas 
fprintf('\n Atividade de Sala 4 \n');
fprintf('          Item a.      \n');
fprintf('\n Alcanse para velocidade de %E [ m/s ]: R = %E [ m ] \n', V, R2.V);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia para velocidade de %E [ m/s ]: E = %E [ s ] \n', V, E2.V);    % [ m/s ] Saída do resultado

fprintf('\n          Item b.      \n');
fprintf('\n Alcanse máximo para velocidade de máximo alcanse %E [ m/s ]: R = %E [ m ] \n', Vmax(1), R2.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia máxima para velocidade de máxima autonomia %E [ m/s ]: E = %E [ s ] \n', Vmax(2), E2.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n ============================================================== \n');    % [ m/s ] Saída do resultado

% ==========================================================================

%% Atividade de sala 📚 @5
% l. Uma aeronave de transporte de passageiros a jato possui as
%     seguintes características:

%                    S = 341,5m2 CD = 0,016 + 0,065CL2; TSFC = 0.8 1/h

%    Para um voo a 30.000 ft, essa aeronave possui massa de 100.000 kg,
%     sendo 30.000 kg de massa de combustível.

%     a) Se a aeronave voa a uma velocidade de 325,8 knot, qual o
%         alcance e a autonomia que ela consegue atingir sem
%         reabastecimento?

%     b) Qual seria a velocidade que a aeronave deveria voar para obter o
%         máximo alcance e a máxima autonomia? E quais seriam os
%         valores do máximo alcance e da máxima autonomia?

% OBS: Considere velocidade e CL constantes.

% Dados 
g = 9.81;                  % [ m/s^2 ] Aceleração da gravidade
S = 341.5;                 % [ m^2 ] Área da asa
m = 100000;                % [ kg ] Massa da aeronave
m_F = 30000;               % [ kg ] Massa de combustível
h = 30000 * 3048;          % [ m ] Altitude
V = 325.8 * 0.514444;      % [ m/s ] Velocidade da aeornave
CD0 = 0.016;              % [ - ] Coeficiente de Arrasto
k = 0.065;                 % [ - ] Coeficiente de Arrasto
TSFC = .8 / 3600;          % [ 1/s ] Consumo específico de combustível 

% a) Alcance e autonomia

% Obtenção de dados de atmosfera
[~, ~, ~, rho] = atmosisa( h );

% Obtenção da velocidade e coeficiente de sustentação de mínimo arrasto
CLmd = sqrt( CD0 / k );
Vmd = sqrt( 2 * m * g  / ( CLmd * rho * S ) );

% Obtenção da Máxima Eficiencia Aerodinâmica
Emax = 1 / ( 2 * sqrt( k * CD0 ) );

% Obtenção da Razão de Velocidade
u = V / Vmd;

% Obtenção da Razão de peso de cruzeiro
ome = ( m * g ) / ( ( m - m_F ) * g );

% Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
R3.V = ( Vmd / TSFC * Emax) * 2 * u * ( atan( 1 / ( u ^2 ) ) - atan( 1 / ( ome * u ^2 ) ) );
[ ~ , index ] = max( R3.V ) ;  

% Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
E3.V = ( 1 / TSFC * Emax) * 2 * ( atan( 1 / ( u ^2 ) ) - atan( 1 / ( ome * u ^2 ) ) );
[ ~ , jndex ] = max( E3.V ) ;  

% Velocidade para máximo alcanse e Máximo Alcanse
Vmax(1) = Vmd * u(index);
u3R.max = Vmax(1) / Vmd;
R3.Vmax = ( Vmd / TSFC * Emax) * 2 * u3R.max * ( atan( 1 / ( u3R.max ^2 ) ) - atan( 1 / ( ome * u3R.max ^2 ) ) );

% Velocidade para máxima autonomia e máxima Autonomia
Vmax(2) = Vmd * u(jndex);
u3E.max = Vmax(2) / Vmd;
E3.Vmax = ( 1 / TSFC * Emax) * 2 * ( atan( 1 / ( u3E.max ^2 ) ) - atan( 1 / ( ome * u3E.max ^2 ) ) );

% Ipressão de respostas 
fprintf('\n Atividade de Sala 4 \n');
fprintf('          Item a.      \n');
fprintf('\n Alcanse para velocidade de %E [ m/s ]: R = %E [ m ] \n', V, R2.V);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia para velocidade de %E [ m/s ]: E = %E [ s ] \n', V, E2.V);    % [ m/s ] Saída do resultado

fprintf('\n          Item b.      \n');
fprintf('\n Alcanse máximo para velocidade de máximo alcanse %E [ m/s ]: R = %E [ m ] \n', Vmax(1), R2.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n Autonomia máxima para velocidade de máxima autonomia %E [ m/s ]: E = %E [ s ] \n', Vmax(2), E2.Vmax);    % [ m/s ] Saída do resultado
fprintf('\n ============================================================== \n');    % [ m/s ] Saída do resultado

% ==========================================================================

%% Atividade de casa 📚 @1
% 1. Considerando omega = 1,5, plote as funções de alcance para os 3 métodos
%     para um intervalor de velocidade relativa de 0,5 a 3,0. Qual método
%     apresenta o melhor alcance? É possível voar nessa condição do ponto
%     de vista de operação?

% Dados 
g = 9.81;                  % [ m/s^2 ] Aceleração da gravidade
S = 341.5;                 % [ m^2   ] Área da asa
m = 100000;                % [ kg    ] Massa da aeronave
m_F = 30000;               % [ kg    ] Massa de combustível
h = 30000 * 3048;          % [ m     ] Altitude
V = 325.8 * .514444;       % [ m/s   ] Velocidade da aeronave 
CD0 = 0.016;               % [ -     ] Coeficiente de Arrasto
k = 0.065;                 % [ -     ] Coeficiente de Arrasto
TSFC = .8 / 3600;          % [ 1/s   ] Consumo específico de combustível 

% a) Alcance e autonomia

% Obtenção de dados de atmosfera
[~, ~, ~, rho] = atmosisa( h );

% Obtenção da velocidade e coeficiente de sustentação de mínimo arrasto
CLmd = sqrt( CD0 / k );
Vmd = sqrt( 2 * m * g  / ( CLmd * rho * S ) );

% Obtenção da Razão de Velocidade e Velocidade
u = linspace(.5, 3, 1E3);

% Obtenção da Máxima Eficiencia Aerodinâmica
Emax = 1 / ( 2 * sqrt( k * CD0 ) );

% Obtenção da Razão de peso de cruzeiro
ome = 1.5;

% MÉTODO 1 ===========================================================================
    % Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
    R1.V = ( Vmd / TSFC * Emax) * ( 2 .* u .^3 ./ ( u .^4 + 1 ) ) * log( ome );

    % Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
    E1.V = ( 1 / TSFC * Emax) .* ( 2 .* u .^2 ./ ( u .^4 + 1 ) ) * log( ome );

    % Velocidade para máximo alcanse e Máximo Alcanse
    Vmax(1) = Vmd * ( 3 ) ^ ( 1 / 4 );
    u1R.max  = Vmax(1) / Vmd;
    R1.Vmax = ( Vmd / TSFC * Emax) * ( 2 * u1R.max ^3 / ( u1R.max ^4 + 1 ) ) * log( ome );

    % Velocidade para máxima autonomia e máxima Autonomia
    Vmax(2) = Vmd * 1;
    u1E.max  = Vmax(2) / Vmd;
    E1.Vmax = ( 1 / TSFC * Emax) * ( 2 * u1E.max ^3 / ( u1E.max ^4 + 1 ) ) * log( ome );
% ===================================================================================

% MÉTODO 2  =========================================================================
    % Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
    R2.V = ( Vmd / TSFC * Emax) * ( 2 .* u .^3 ./ ( u .^4 + 1 ) ) * 2 * ( 1 - sqrt( ome ) ^( -1 ) );

    % Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
    E2.V = ( 1 / TSFC * Emax) * ( 2 * u .^2 ./ ( u .^4 + 1 ) ) * log( ome );

    % Velocidade para máximo alcanse e Máximo Alcanse
    Vmax(1) = Vmd * ( 3 ) ^ ( 1 / 4 );
    u2R.max = Vmax(1) / Vmd;
    R2.Vmax = ( Vmd / TSFC * Emax) * ( 2 * u2R.max ^3 / ( u2R.max ^4 + 1 ) ) * 2 * ( 1 - sqrt( ome ) ^( -1 ) );

    % Velocidade para máxima autonomia e máxima Autonomia
    Vmax(2) = Vmd * 1;
    u2E.max  = Vmax(2) / Vmd;
    E2.Vmax = ( 1 / TSFC * Emax) * ( 2 * u2E.max ^2 / ( u2E.max ^4 + 1 ) ) * log( ome );
% ==========================================================================


% MÉTODO 3  =========================================================================
    % Obtenção do Alcanse para velocidade V = 325.8 [ knot ]
    R3.V = ( Vmd / TSFC * Emax) * 2 .* u .* ( atan( 1 ./ ( u .^2 ) ) - atan( 1 ./ ( ome .* u .^2 ) ) );
    [ ~ , index ] = max( R3.V ) ;  

    % Obtenção da Autonomia para velocidade V = 325.8 [ knot ]
    E3.V = ( 1 / TSFC * Emax) * 2 .* ( atan( 1 ./ ( u .^2 ) ) - atan( 1 ./ ( ome .* u .^2 ) ) );
    [ ~ , jndex ] = max( E3.V ) ;  

    % Velocidade para máximo alcanse e Máximo Alcanse
    Vmax(1) = Vmd * u(index);
    u3R.max = Vmax(1) / Vmd;
    R3.Vmax = ( Vmd / TSFC * Emax) * 2 * u3R.max * ( atan( 1 / ( u3R.max ^2 ) ) - atan( 1 / ( ome * u3R.max ^2 ) ) );

    % Velocidade para máxima autonomia e máxima Autonomia
    Vmax(2) = Vmd * u(jndex);
    u3E.max = Vmax(2) / Vmd;
    E3.Vmax = ( 1 / TSFC * Emax) * 2 * ( atan( 1 / ( u3E.max ^2 ) ) - atan( 1 / ( ome * u3E.max ^2 ) ) );
% ==========================================================================

% Gráfico de Alcance
figure;hold on; grid minor;
plot( u, R1.V / 1000, 'r','LineWidth', 1 )
plot( u, R2.V / 1000, 'k','LineWidth', 1 )
plot( u, R3.V / 1000, 'b','LineWidth', 1 )
plot( u1R.max, R1.Vmax / 1000, 'ro','LineWidth', 2 )
xline( u1R.max, 'k-.', 'u1', 'LineWidth', .5 )
plot( u2R.max, R2.Vmax / 1000, 'ko','LineWidth', 2 )
xline( u2R.max, 'k-.', 'u2', 'LineWidth', .5 )
plot( u3R.max, R3.Vmax / 1000, 'bo','LineWidth', 2 )
xline( u3R.max, 'k-.', 'u3', 'LineWidth', .5 )
title('Alcance','FontSize', 12,'FontName','Times New Roman');
xlabel('Velocidade Relativa [ - ]','FontSize', 12,'FontName','Times New Roman');
ylabel('Alcance [ km ]','FontSize', 12,'FontName','Times New Roman');
legend('Método 1', 'Método 2', 'Método 3', 'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')
% ===========================================================================================================

%% Atividade de casa 📚 @2
% 2. Considerando omega = 1,5, plote as funções de autonomia para os 3
%     métodos para um intervalor de velocidade relativa de 0,5 a 3,0. Qual
%     método apresenta o melhor autonomia? É possível voar nessa condição
%     do ponto de vista de operação?

% Gráfico de Autonomia
figure;hold on; grid minor;
plot( u, E1.V / 3600, 'r', 'LineWidth', 2 )
plot( u, E2.V / 3600, 'k--', 'LineWidth', 1 )
plot( u, E3.V / 3600, 'b', 'LineWidth', 1 )
plot( u1E.max, E1.Vmax / 3600, 'ro','LineWidth', 2 )
xline( u1E.max, 'k-.', 'u1', 'LineWidth', .5 )
plot( u2E.max, E2.Vmax / 3600, 'ko','LineWidth', 2 )
xline( u2E.max, 'k-.', 'u2', 'LineWidth', .5 )
plot( u3E.max, E2.Vmax / 3600, 'bo','LineWidth', 2 )
xline( u3E.max, 'k-.', 'u3', 'LineWidth', .5 )
title('Autonomia','FontSize', 12,'FontName','Times New Roman');
xlabel('Velocidade Relativa [ - ]','FontSize', 12,'FontName','Times New Roman');
ylabel('Autonomia [ h ]','FontSize', 12,'FontName','Times New Roman');
legend('Método 1', 'Método 2', 'Método 2', 'FontSize', 12,'FontName','Times New Roman','Location','best');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); set(gcf,'paperPositionMode','auto')

% ===========================================================================================================

%% Atividade de casa 📚 @3
% 3. Faça um resumo das reservas técnicas de combustível definidas pelo
%     RBAC NO 121.645 após a Emenda no 16.


