%% Atividade 1 - Desenpenho
% Professor: Higor Luis
clc; close all; clear all; format compact
% ========================================================================

% ========================================================================
%%
fprintf("\n ============================================== \n")
fprintf("\n                 QUESTÃO 1                      \n")
fprintf("\n ============================================== \n")

% 1. Plotar as propriedades relativas para uma variação de 
%    altitude de 0 a 100.000 ft.

fprintf("\n ________________________________________ \n")
fprintf("\n Impressão de Gráficos \n")
fprintf("\n ________________________________________ \n")

% Variáveis
h1   = 100000 * .3048;          % Transformação ft -> m 
H1   = linspace(0, h1, 200000); % Altitude variando [m]

% Valores Padrão
rho = 1.225;
T   = 288.15;
P   = 101325;

% Função atmosisa
[T1, a1, P1, rho1] = atmosisa(H1);

% Gráficos
figure
subplot(2,2,1)
plot(rho1,H1);title('Densidade vs. Altitude','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
xlabel('Densidade $[\frac{kg}{m^3}]$','Interpreter','latex','FontSize', 12,'FontName','Times New Roman'); ylabel('Altitude [m]','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
legend('Densidade','FontSize', 12,'FontName','Times New Roman',...
       'Location','best','Interpreter','latex');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
set(gcf,'paperPositionMode','auto');

subplot(2,2,2)
plot(T1,H1);title('Temperatura vs. Altitude','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
xlabel('Temperatura $[K]$','Interpreter','latex','FontSize', 12,'FontName','Times New Roman'); ylabel('Altitude [m]','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
legend('Densidade','FontSize', 12,'FontName','Times New Roman',...
       'Location','best','Interpreter','latex');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
set(gcf,'paperPositionMode','auto');

subplot(2,2,3)
plot(P1,H1);title('Pressão vs. Altitude','FontSize', 12,'FontName','Times New Roman')
xlabel('Pressão [Pa]','FontSize', 12,'FontName','Times New Roman'); ylabel('Altitude [m]','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
legend('Densidade','FontSize', 12,'FontName','Times New Roman',...
       'Location','best','Interpreter','latex');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
set(gcf,'paperPositionMode','auto');

subplot(2,2,4)
plot(a1,H1);title('Velocidade do Som vs. Altitude','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
xlabel('Velocidade do som $[\frac{m}{s}]$','Interpreter','latex','FontSize', 12,'FontName','Times New Roman'); ylabel('Altitude [m]','Interpreter','latex','FontSize', 12,'FontName','Times New Roman')
legend('Densidade','FontSize', 12,'FontName','Times New Roman',...
       'Location','best','Interpreter','latex');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
set(gcf,'paperPositionMode','auto')

figure 
hold on
plot (rho1/rho,H1);
plot (T1/T,H1);
plot (P1/P,H1);
xlabel('Propriedades Relativas','FontSize', 12,'FontName','Times New Roman'); 
ylabel('Altitude [m]','FontSize', 12,'FontName','Times New Roman')
legend('$\sigma$','$\theta$','$\delta$','Interpreter','latex','FontSize', 12,'FontName','Times New Roman',...
       'Location','best','Interpreter','latex');
set(gcf, 'Color', 'w'); set(gca,'GridLineStyle', '-'); 
set(gcf,'paperPositionMode','auto')


% ========================================================================
%%
fprintf("\n ============================================== \n")
fprintf("\n                 QUESTÃO 2                      \n")
fprintf("\n ============================================== \n")

% 2. Determine a temperatura, pressão, densidade e velocidade 
%    do som a 5000 m e condição ISA. Recalcule para 
%    a mesma altitude e ISA-10.

% Variáveis
h2   = 5000;                     % Altitude [m] 
H2  = linspace(0, h2*2, 200000); % Altitude variando [m]

% Função atmosisa
[T2, a2, P2, rho2] = atmosisa(h2);
[T2var, a2var, P2var, rho2var] = atmosisa(H2);

% Impressão de dados ISA para H2 = 5000 m
fprintf("\n ________________________________________ \n")
fprintf("\n H2 = 5000 [m] \n")
    fprintf("\n Temperatura = %G [K] \n", T2)
    fprintf("\n Velocidade do som = %G [m/s] \n", a2)
    fprintf("\n Pressão = %G [Pa] \n", P2)
    fprintf("\n Densidade = %G [kg/m^3] \n", rho2)
fprintf("\n ________________________________________ \n")
    
% Obtenção de dados para ISA-10
fprintf("\n ________________________________________ \n")
fprintf("\n A temperatura deveria ser %G°C e o obtido foi %G°C \n", (T2-273.15), (T2+10)-273.15)

% Localização dos dados
[~, index] = min( abs( T2var - ( T2 + 10 ) ) );
rho2       = rho2var(index);
P2         = P2var(index);
a2         = a2var(index);

% Impressão de dados ISA-10
fprintf("\n Dados para ISA-10 \n")
    fprintf("\n A Temperatura para ISA-10 = %G [K] \n", T2+10)
    fprintf("\n A Densidade para ISA-10 =  %G [kg/m^3] \n", rho2)
    fprintf("\n Velocidade do som para ISA-10 = %G [m/s] \n", a2)
    fprintf("\n Pressão para ISA-10 = %G [Pa] \n", P2)
fprintf("\n ________________________________________ \n")


% ========================================================================
%%
fprintf("\n ============================================== \n")
fprintf("\n                 QUESTÃO 3                      \n")
fprintf("\n ============================================== \n")

% 3. Uma aeronave voa a uma altitude onde a temperatura é -4,5 °C. 
%    Calcule a altitude na condição ISA e ISA+10.

% Variáveis
h3  = 100000;              % Altitude [m] 
H3  = linspace(0,h3,200000); % Altitude variando [m]

% Função atmosisa
[T3var, a3var, P3var, rho3var] = atmosisa(H3);

% Localização dos dados
[~, index] = min( abs( ( -4.5 + 273.15 ) - ( T3var ) ) );
H3isa      = H3(index);
rho3       = rho3var(index);
P3         = P3var(index);
a3         = a3var(index);
T3         = T3var(index);

% Impressão de dados para T = -4.5 °C
fprintf("\n ________________________________________ \n")
fprintf("\n Impressão de dados para T = -4.5 °C \n")
    fprintf("\n Altura = %G [m] \n", H3isa)
    fprintf("\n Temperatura = %G [K] \n", T3)
    fprintf("\n Velocidade do som = %G [m/s] \n", a3)
    fprintf("\n Pressão = %G [Pa] \n", P3)
    fprintf("\n Densidade = %G [kg/m^3] \n", rho3)
fprintf("\n ________________________________________ \n")

    % Localização dos dados
[~, index] = min( abs( ( T3var ) - ( -4.5 - 10 + 273.15 ) ) );
rho3       = rho3var(index);
P3         = P3var(index);
a3         = a3var(index);
T3         = T3var(index);
H3         = H3(index);

% Impressão de dados para T = 5.5 °C
fprintf("\n ________________________________________ \n")
fprintf("\n Impressão de dados para T = -14.5 °C \n")
    fprintf("\n Altura = %G [m] \n", H3)
    fprintf("\n Temperatura = %G [K] \n", T3)
    fprintf("\n Velocidade do som = %G [m/s] \n", a3)
    fprintf("\n Pressão = %G [Pa] \n", P3)
    fprintf("\n Densidade = %G [kg/m^3] \n", rho3)
fprintf("\n ________________________________________ \n")


% ========================================================================
%%

fprintf("\n ============================================== \n")
fprintf("\n                 QUESTÃO 3                      \n")
fprintf("\n ============================================== \n")
% 4. Qual é a pressão dinâmica quando uma aeronave voa em 
%    cruzeiro a uma altitude de 14500 m e Mach 0,6 e ISA-12?


% Variáveis
h4  = 14500;              % Altitude [m] 
H4  = linspace(0,h4*2,20000); % Altitude variando [m]

% Função atmosisa
[T4, a4, P4, rho4] = atmosisa(h4);
[T4var, a4var, P4var, rho4var] = atmosisa(H4);

% Impressão de dados ISA para H2 = 14500 m
fprintf("\n ________________________________________ \n")
fprintf("\n H2 = 14500 [m] \n")
    fprintf("\n Temperatura = %G [K] \n", T4)
    fprintf("\n Velocidade do som = %G [m/s] \n", a4)
    fprintf("\n Pressão = %G [Pa] \n", P4)
    fprintf("\n Densidade = %G [kg/m^3] \n", rho4)
    fprintf("\n Pressão dinamica = %G [N/m^2] \n", .5*rho4*(a4*0.6)^2)
fprintf("\n ________________________________________ \n")

fprintf("\n ________________________________________ \n")
% Obtenção de dados para ISA-12
fprintf("\n A temperatura deveria ser %G°C, mas a presente era de %G°C \n", (T4-273.15), (T4+12)-273.15)

% Localização dos dados
[~, index] = min( abs( ( T4var ) - ( T4 + 12  ) ) );
rho4       = rho4var(index);
P4         = P4var(index);
a4         = a4var(index);
T4         = T4var(index);


% Impressão de dados para ISA-12
fprintf("\n Impressão de dados para ISA-12 \n")
    fprintf("\n Temperatura = %G [K] \n", T4)
    fprintf("\n Velocidade do som = %G [m/s] \n", a4)
    fprintf("\n Pressão = %G [Pa] \n", P4)
    fprintf("\n Densidade = %G [kg/m^3] \n", rho4)
    fprintf("\n Pressão dinamica = %G [N/m^2] \n", .5*rho4*(a4*0.6)^2)
fprintf("\n ________________________________________ \n")

    
    
