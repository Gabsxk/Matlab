%% Inicialização 🤔🛠️
clc; close all; clear all; format compact;
% =======================================================================

%% Atividade de Sala 1
% A aeronave Beech Baron 58 está voando ao nível do mar com uma
%  velocidade de 80 m/s e 3° de ângulo de ataque. A área da asa é de
%  18,51 m2 e seu alongamento (AR) é de 7,2. Determine as forças
%  aerodinâmicas (sustentação e arrasto) que esta aeronave está
%  produzindo. Assuma que o coeficiente de arrasto da aeronave é 0,05 e o
%  ângulo de ataque de sustentação zero da asa é zero (a0 = 0).

% Dados
V = 80;        % [ m/s ] Velocidade
h = 0;         % [ m   ] Altitude
S = 18.51;     % [ m2  ] Área da asa
AR = 7.2;      % [  -  ] Alongamento
aw = 3*pi/180; % [ rad ] Ângulo de ataque
a0 = 2*pi;     % [ rad ] Aplicando a teoria dos aerofólios finos
af = 0;        % [ rad ] Ângulo de ataque de sustentação zero
CD0 = 0.05;    % [  -  ] Coeficiente de arrasto

% Cálculos
[ T, a, P, rho ] = atmosisa( h );     % Obtenção dos dados atmosféricos
CLaw = a0 / ( 1 + ( a0 / (pi*AR) ) ); % Inclinação da curva de sustentação
CLw  = CLaw *  (aw - af);             % Coeficiente de sustentação
L = 0.5 * rho * V^(2) * S * CLw;      % [ N ] Força de sustentação
D = 0.5 * rho * V^(2) * S * CD0;      % [ N ] Força de arrasto

% Impressão
fprintf( '\n Atividade 1 de sala \n')
fprintf( '\n Força de sustentação: %E \n', L) % [ N ]
fprintf( '\n Força de arrasto: %E \n', D)     % [ N ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Atividade de Sala 2
% Uma aeronave com massa de 500 kg e área de asa de 8 m2 está voando
%  ao nível do mar com velocidade constante de 50 m/s. O motor da
%  aeronave está gerando 600 N de empuxo. Determine o coeficiente de
%  sustentação e o coeficiente de arrasto da aeronave nesta condição de
%  voo. Assuma um ângulo de ataque nulo. Em seguida, determine a razão
%  sustentação-arrasto.

% Dados
m  = 500;       % [ kg  ] Massa
V  = 50;        % [ m/s ] Velocidade
h  = 0;         % [ m   ] Altitude
S  = 8;         % [ m2  ] Área da asa
Tr = 600;       % [ N   ] Empuxo
aw = 0;         % [ rad ] Ângulo de ataque

% Cálculos
[ T, a, P, rho ] = atmosisa( h );                % Obtenção dos dados atmosféricos
CLw  = ( m * 9.81 ) / ( 0.5 * rho * V^(2) * S ); % Coeficiente de sustentação
L = 0.5 * rho * V^(2) * S * CLw;                 % [ N ] Força de sustentação
D = Tr * cos(a);                                % [ N ] Força de arrasto
CD = D / ( 0.5 * rho * V^(2) * S );              % [  - ] Coeficiente de arrasto
LD = L / D;                                      % [  - ] Razão sustentação-arrasto

% Impressão
fprintf( '\n Atividade 2 de sala \n')
fprintf( '\n Coeficiente de sustentação: %E \n', CLw) % [  - ]
fprintf( '\n Coeficiente de arrasto: %E \n', CD)      % [  - ]
fprintf( '\n Razão sustentação-arrasto: %E \n', LD)   % [  - ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Atividade de Sala 3
% A aeronave da figura está subindo com um ângulo de subida de 15° e
%  um ângulo de ataque de 5°. Se a aeronave tem uma massa de 12.000 kg
%  e o arrasto da aeronave é de 10.000 N, determine o empuxo do motor e
%  a sustentação da aeronave.

% Dados
m   = 12000;     % [ kg  ] Massa
V   = 50;        % [ m/s ] Velocidade
D   = 10000;     % [ N   ] Arrasto
gam = 15*pi/180; % [ rad ] Ângulo de subida
aw  = 5*pi/180;  % [ rad ] Ângulo de ataque da fuselagem

% Cálculos
Tr = ( D + ( m * 9.81 ) * sin( gam ) ) / cos(aw); % [ N ] Empuxo
L  = ( m * 9.81 ) * cos( gam ) - Tr * sin( aw );  % [ N ] Força de sustentação

% Impressão
fprintf( '\n Atividade 3 de sala \n')
fprintf( '\n Empuxo: %E N\n', Tr) % [ N ]
fprintf( '\n Sustentação: %E N\n', L) % [ N ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Atividade de Sala 4
% O caça F-15 tem uma massa de 30.845 kg e uma área de asa de 56,5 m2.
%  Se este caça estiver voando a 15.000 m de altitude com o coeficiente de
%  sustentação de 0,1, determine sua velocidade verdadeira e equivalente
%  em knots

% Dados
hp  = 0;         % [ m   ] Altitude Padrão
m   = 30845;     % [ kg  ] Massa
S   = 56.5;      % [ m2  ] Área da asa
h   = 15000;     % [ m   ] Altitude
CL  = 0.1;       % [  -  ] Coeficiente de sustentação

% Cálculos
[ T, a, P, rho ] = atmosisa( [ h hp] );           % Obtenção dos dados atmosféricos
Vt = sqrt( ( 2 * m * 9.81 ) / ( rho(1) * S * CL ) ); % [ m/s ] Velocidade verdadeira
Ve = Vt / sqrt( rho(1) / rho(2)  );               % [ m/s ] Velocidade equivalente

% Impressão
fprintf( '\n Atividade 4 de sala \n')
fprintf( '\n Velocidade verdadeira: %E m/s\n', Vt) % [ m/s ]
fprintf( '\n Velocidade equivalente: %E m/s\n', Ve) % [ m/s ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Atividade de Sala 5
% Considere uma aeronave de carga com as seguintes caraterísticas:
%  m = 380.000 kg; S = 576 m2; MAC = 9.3;  tcm = .18; Cdmin = .0052
%  Esta aeronave está voando ao nível do mar com uma velocidade de 400
%  knot. Assumindo que o CD o da aeronave é três vezes o CDO da asa (isto é,
%  CDow), determine o CDO da aeronave.

% Dados
m     = 380000;   % [ kg  ] Massa
S     = 576;      % [ m2  ] Área da asa
MAC   = 9.3;      % [ m   ] Comprimento médio da corda
tcm   = .18;      % [  -  ] Espessura relativa da corda média
Cdmin = .0052;    % [  -  ] Coeficiente de arrasto mínimo
V     = 400 * 0.514444; % [ m/s ] Velocidade
h     = 0;        % [ m    ] Altitude
mu    = 1.81e-5;  % [ kg/ms] Viscosidade dinâmica do ar

% Cálculos
[ T, a, P, rho ] = atmosisa( h ); % Obtenção dos dados atmosféricos
Re = rho * V * MAC / mu; % [  - ] Número de Reynolds
M  = V / a;              % [  - ] Número de Mach

% Análise de Regime de Escoamento
fprintf( '\n Atividade 5 de sala \n')
if Re < 2e5
    fprintf( '\n Regime laminar \n')
else
    fprintf( '\n Regime turbulento \n')
end

% Cálculo do CDO
Cf   = .455 / ( log10( Re ) )^( 2.58 );
fm   = 1 - .08 * M ^ ( 1.45 );
ftc  = 1 + 2.7 * tcm + 100 * tcm ^ ( 4 );
Swet = 2 * ( 1 + 0.5 * tcm ) * S;
Cd0w = Cf * fm * ftc * Swet / S * ( Cdmin / .004 ) ^ ( 0.4 );
Cd0  = 3 * Cd0w ; 

% Impressão
fprintf( '\n Coeficiente de arrasto da aeornave: %E \n', Cd0)
fprintf( '===============================================')
% =======================================================================

%% Exercício de Entrega 1
% 1.Uma aeronave de carga com peso de 145.000 lb e área de asa de 1.318 ft2
% tem um coeficiente de sustenta ao máximo de 2,5. Esta aeronave é capaz
% de voar a uma altitude de 55.000 ft e condição ISA+15 com uma
% velocidade de 150 KTAS?

% Dados
hp  = linspace( 0, 80000, 1000000); % [ m   ] Altitude Padrão
m   = 145000 * .453592; % [ kg  ] Massa
S   = 1318 * .092903; % [ m2  ] Área da asa
h   = 55000 * .3048; % [ m   ] Altitude
CL  = 2.5; % [  -  ] Coeficiente de sustentação
V   = 150 * .514444; % [ m/s ] Velocidade verdadeira

% Cálculos
[ T, a, P, rho ] = atmosisa( h ); % Obtenção dos dados atmosféricos

% Análise da Condição ISA+15
[ Tp, ap, Pp, rhop ] = atmosisa( hp ); % Obtenção dos dados atmosféricos
fprintf( '\n Atividade 1 de Entrega \n')
fprintf( '\n A temperatura deveria ser %E°C e o obtido foi %E°C \n', ( T - 273.15 ), ( T - 15 ) - 273.15)

% Localização dos dados
[~, index] = min( abs( Tp - ( T - 15 ) ) );
rho        = rhop(index);
P          = Pp( index );
a          = ap( index );

% Análise de possibilidade de voo
L = 0.5 * rho * V^2 * S * CL; % [ N ] Força de sustentação

if L > m * 9.81
    fprintf( '\n A aeronave pode voar L > W \n')
else
    fprintf( '\n A aeronave não pode voar L < W\n')
end
fprintf( '-----------------------------------------')
% =======================================================================

%% Exercício de Entrega 2
% 2. O bombardeiro B-IB tem massa máxima de decolagem de 216.367 kg,
% área de asa de 181 m2 e velocidade máxima de Mach 2,2. Supondo que o
% coeficiente de arrasto desta aeronave em cruzeiro seja 0,03, quanto
% empuxo os quatro motores geram para esta condição de voo?

% Dados
h   = 0;                            % [ m   ] Altitude
m   = 216367;                       % [ kg  ] Massa
S   = 181;                          % [ m2  ] Área da asa
M   = 2.2;                          % [  -  ] Mach
CD  = 0.03;                         % [  -  ] Coeficiente de arrasto

% Cálculos
[ T, a, P, rho ] = atmosisa( h );              % Obtenção dos dados atmosféricos
V = M * a;                                     % [ m/s ] Velocidade verdadeira
Tr = 0.5 * rho * V^2 * S * CD;                 % [ N ] Empuxo

% Impressão
fprintf( '\n Atividade 2 de Entrega \n')
fprintf( '\n Empuxo: %E N \n', 4*Tr) % [ N ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Exercício de Entrega 3
% 3. A aeronave Voyager é capaz de voar ao redor do globo sem reabastecer.
% Em uma missão, a aeronave está voando no equador a uma altitude de
% 15.000 pés com uma velocidade de 110 knot. Suponha que haja um vento
% de 15 m/s soprando de oeste para leste o tempo todo.
% Nota: A Terra tem um diâmetro de 12.800 km.

% Dados
h   = 15000*0.3048;                 % [ m   ] Altitude
V   = 110*0.514444;                 % [ m/s ] Velocidade verdadeira
Vw  = 15;                           % [ m/s ] Velocidade do vento
Dm  = 12800*1000;                   % [ m   ] Diâmetro da Terra

% a.Quantos dias leva para fazer esta missão se estiver navegando de oeste
% para leste?

% Cálculos
[ T, a, P, rho ] = atmosisa( h );   % Obtenção dos dados atmosféricos
Vt = V + Vw;                        % [ m/s ] Velocidade verdadeira
t  = Dm / Vt;                       % [ s   ] Tempo de voo

% Impressão
fprintf( '\n Atividade 3 de Entrega \n')
fprintf( '\n a. Voando de oeste para leste \n')
fprintf( '\n %E dias \n', t/86400) % [ dia ]

% b.Quantos dias leva para fazer esta missão se estiver navegando do leste
% para o oeste?

% Cálculos
[ T, a, P, rho ] = atmosisa( h );   % Obtenção dos dados atmosféricos
Vt = V - Vw;                        % [ m/s ] Velocidade verdadeira
t  = Dm / Vt;                       % [ s   ] Tempo de voo

fprintf( '\n b. Voando de leste para oeste \n')
fprintf( '\n %E dias \n', t/86400) % [ dia ]
fprintf( '-----------------------------------------')
% =======================================================================

%% Exercício de Entrega 4
% 4.Uma aeronave GA está voando a 5.000 pés de altitude com uma
%    velocidade de 100 knot. O comprimento da fuselagem é de 7 m, o
%    MAC da asa é de 1,5 m, o MAC da cauda horizontal é de 0,8 m e o
%    MAC da cauda vertical é de 0,6 m. Determine o número de Reynolds
%    da fuselagem, asa, empenagens horizontal e vertical.

% Dados
h   = 5000*0.3048;                  % [ m   ] Altitude
V   = 100*0.514444;                 % [ m/s ] Velocidade verdadeira
L   = 7;                            % [ m   ] Comprimento da fuselagem
MACa = 1.5;                         % [ m   ] MAC da asa
MACh = 0.8;                         % [ m   ] MAC da cauda horizontal
MACv = 0.6;                         % [ m   ] MAC da cauda vertical

% Obtenção dos dados atmosféricos
[ T, a, P, rho ] = atmosisa( h );

% Determinando viscosidade dinamica
mu = 1.458e-6 * ( T ^ ( 1 / 2 ) ) / ( 1 + 110.4 / T ); % [ kg/(m.s) ] Viscosidade dinâmica Eq. 1.26 p 18

% Determinando número de Reynolds
Re_f = rho * V * L / mu;            % [  -  ] Número de Reynolds da fuselagem
Re_a = rho * V * MACa / mu;         % [  -  ] Número de Reynolds da asa
Re_h = rho * V * MACh / mu;         % [  -  ] Número de Reynolds da cauda horizontal
Re_v = rho * V * MACv / mu;         % [  -  ] Número de Reynolds da cauda vertical

% Impressão
fprintf( '\n Atividade 4 de Entrega \n')
fprintf( '\n Número de Reynolds da fuselagem: %E \n', Re_f)
fprintf( '\n Número de Reynolds da asa: %E \n', Re_a)
fprintf( '\n Número de Reynolds da cauda horizontal: %E \n', Re_h)
fprintf( '\n Número de Reynolds da cauda vertical: %E \n', Re_v)
fprintf( '\n ----------------------------------------- \n')
% =======================================================================

%% Exercício de Entrega 5
% 5. A aeronave Moony M20TN (Figura no próximo slide) possui uma asa
%     com as seguintes características:
%             S = 16,3 m2, b = 11,11 m, (t/c)max = 15%,
%               Airfoil: NACA 63-215; Cdminw = 0,0042 
%    A aeronave tem uma massa de 1528 kg e está voando a uma
%     velocidade de 237 knot. 

% Dados
S     = 16.3;         % [ m2  ] Área da asa
b      = 11.11;        % [ m   ] Envergadura
tcm    = 0.15;         % [  -  ] Espessura relativa máxima
m      = 1528;         % [ kg  ] Massa da aeronave
V      = 237*0.514444; % [ m/s ] Velocidade verdadeira
L      = 8.128;        % [ m   ] Comprimento da fuselagem
Cdminw = 0.0042;       % [  -  ] Coeficiente de arrasto mínimo
h      = 0;            % [ m   ] Altitude - Assumido

% Obtenção das cordas médias aerodinâmicas
MACa   = ( 2 / 3 ) * 2.0 * (1 + ( .10 / 2.0 ) - ( ( .10 / 2.0 ) / ( 1 + ( .10 / 2.0 ) ) ) ); % [ m   ] MAC da asa
MACh   = ( 2 / 3 ) * .93 * (1 + ( .62 / .93 ) - ( ( .62 / .93 ) / ( 1 + ( .62 / .93 ) ) ) ); % [ m   ] MAC da cauda horizontal
MACv   = ( 2 / 3 ) * .12 * (1 + ( .62 / .12 ) - ( ( .62 / .12 ) / ( 1 + ( .93 / .62 ) ) ) ); % [ m   ] MAC da cauda vertical

% Obtenção dos dados atmosféricos
[ T, a, P, rho ] = atmosisa( h );   % Obtenção dos dados atmosféricos

% Determinando viscosidade dinamica
mu = 1.458e-6 * ( T ^ ( 1 / 2 ) ) / ( 1 + 110.4 / T ); % [ kg/(m.s) ] Viscosidade dinâmica Eq. 1.26 p 18

% Determinando número de Reynolds
Re_f = rho * V * L / mu;            % [  -  ] Número de Reynolds da fuselagem
Re_a = rho * V * MACa / mu;         % [  -  ] Número de Reynolds da asa
Re_h = rho * V * MACh / mu;         % [  -  ] Número de Reynolds da cauda horizontal
Re_v = rho * V * MACv / mu;         % [  -  ] Número de Reynolds da cauda vertical

% Análise de Regime de Escoamento
fprintf( '\n Atividade 5 de Entrega \n')
if Re_f & Re_a & Re_h & Re_v < 2e5
    fprintf( '\n Regime Laminar \n')
else
    fprintf( '\n Regime Turbulento \n')
end

% Cálculo do CDOw
M     = V / a;  % [  - ] Número de Mach
Cf    = .455 / ( log10( Re_a ) )^( 2.58 ); % [  -  ] Coeficiente de arrasto da fuselagem
fm    = 1 - .08 * M ^ ( 1.45 ); % [  -  ] Fator de compressibilidade da fuselagem
ftc   = 1 + 2.7 * tcm + 100 * tcm ^ ( 4 ); % [  -  ] Fator de espessura relativa da asa
Swetw = 2 * ( 1 + 0.5 * tcm ) * S; % [ m2  ] Área de superfície molhada da asa
Cd0w  = Cf * fm * ftc * Swetw / S * ( Cdmin / .004 ) ^ ( 0.4 ); % [  -  ] Coeficiente de arrasto de sustentação zero da asa

% Cálculo do CDOh
Sh    = .006392; % [ m2  ] Área da cauda horizontal 
Cf    = .455 / ( log10( Re_h ) )^( 2.58 ); % [  -  ] Coeficiente de arrasto da fuselagem
fm    = 1 - .08 * M ^ ( 1.45 ); % [  -  ] Fator de compressibilidade da fuselagem
ftc   = 1 + 2.7 * tcm + 100 * tcm ^ ( 4 ); % [  -  ] Fator de espessura relativa da cauda horizontal
Sweth = 2 * ( 1 + 0.5 * tcm ) * Sh; % [ m2  ] Área de superfície molhada da cauda horizontal
Cd0h  = Cf * fm * ftc * Sweth / S * ( Cdmin / .004 ) ^ ( 0.4 ); % [  -  ] Coeficiente de arrasto de sustentação zero da asa

% Cálculo do CDOv
Sv    = .006606; % [ m2  ] Área da cauda vertical
Cf    = .455 / ( log10( Re_h ) )^( 2.58 ); % [  -  ] Coeficiente de arrasto da fuselagem
fm    = 1 - .08 * M ^ ( 1.45 ); % [  -  ] Fator de compressibilidade da fuselagem
ftc   = 1 + 2.7 * tcm + 100 * tcm ^ ( 4 ); % [  -  ] Fator de espessura relativa da cauda vertical
Swetv = 2 * ( 1 + 0.5 * tcm ) * Sv; % [ m2  ] Área de superfície molhada da cauda vertical
Cd0v  = Cf * fm * ftc * Swetv / S * ( Cdmin / .004 ) ^ ( 0.4 ); % [  -  ] Coeficiente de arrasto de sustentação zero da asa

% Cálculo do Cdof
fM    = 1 - .08 * M ^ 1.45;
LD    = 812.8 * 1E-2 / 1.8; % [ - ] Razão entre comprimento da aeronave e diâmetro da fuselagem
fLD   = 1 + ( 60 / ( LD ^ 3 ) ) + .0025 * LD; % [  -  ] Fator de fuselagem
Swetf = pi * ( 1.8 / 2 ) ^ ( 2 ) * ( 812.8 * 1E-2 ); % [ m2  ] Área de superfície molhada da asa
Cf    = .455 / ( log10( Re_f ) ) ^ ( 2.58 ); % [  -  ] Coeficiente de arrasto da fuselagem
Cdof  = fM * fLD * Cf * Swetf / S; % [  -  ] Coeficiente de arrasto de sustentação zero da fuselagem

% (1) Determine o coeficiente de arrasto de sustentação zero da asa.

fprintf( '\n (1). Coeficiente de arrasto de sustentação zero da asa \n')
fprintf( '\n Cd0w = %E \n', Cd0w)

% (2) Determine o coeficiente de arrasto de sustentação
%     zero da fuselagem. Para outras informações, use as três 
%     visualizações da aeronave fornecidas.

fprintf( '\n (2). Coeficiente de arrasto de sustentação zero da fuselagem \n')
fprintf( '\n Cd0f = %E \n', Cdof ) % [ dia ]
fprintf('\n ====================== FIM %s \n',char( double( '👍😎' ) ) )
% =======================================================================