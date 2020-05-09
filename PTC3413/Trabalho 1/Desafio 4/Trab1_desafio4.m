format short;
close all;
clc;
global A22 B2 C2
%% Definição dos parâmetros que serão utilizados
t = 0:.1:10; % Tempo de simulação 
M = 4; K = 9; %sqrt(K/M) = 1.5
Tf = 5; % Tempo no qual a entrada age

%% Condições iniciais e finais dos estados de x
x10 = 3; x20 = 1;
x1f = 7; x2f = 7;

%% Definição das matrizes do espaço de estados e das condições iniciais e finais
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,2);

x0 = [x10; 0; x20; 2];
xf = [x1f; 0; x2f; 0];

%% Mudando a base do problema: x = Hz => z = H_inv*x
H_inv = [0.5 0 0.5 0; 0 0.5 0 0.5; 0.5 0 -0.5 0; 0 0.5 0 -0.5];
H = inv(H_inv);

A_til = H_inv*A*H;
B_til = H_inv*B;
C_til = [0.5 -0.5]*(C*H);

z0 = H_inv*x0;
zf = H_inv*xf;

%% Realização da decomposição em partes obsv. e ñ obsv.
Obs = obsv(A_til,C_til);
k = rank(Obs);
[U,S,V] = svd(Obs);
T = [V(:,k+1:end) V(:,1:k)]; % Transformação linear para obter a decomposição

Ad = inv(T)*A_til*T;
Bd = inv(T)*B_til;
Cd = C_til*T;

A11 = Ad(1:k,1:k);
A12 = Ad(1:k,k+1:end);
A22 = Ad(k+1:end,k+1:end);

B1 = Bd(1:k,:);
B2 = Bd(k+1:end,:);

C2 = Cd(:,k+1:end);

z0_obs = z0(k+1:end,1);
z0_uno = z0(1:k,1);

%% Cálculo do grammiano de obsv.
invWo = inv(obsv_gramm(A22,C2,Tf));

%% Simulação do sistema de controle
sim('Trab1_desafio_4_f',10);

%% Demonstração dos resultados

figure('Name', 'Estados estimados pelo sistema de controle');
subplot(2,2,1)
plot(t_sim, z2_calc, 'r', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Estado $z_{2}$ estimado ao longo do tempo', 'interpreter', 'latex');
grid minor;

subplot(2,2,2)
plot(t_sim, erro_z2, 'r--', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Erro do estado $z_{2}$ estimado ao longo do tempo', 'interpreter', 'latex');
grid minor;

subplot(2,2,3)
plot(t_sim, z2p_calc, 'k', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro', 'Interpreter', 'latex');
title('Valor do estado $\dot{z_{2}}$ estimado ao longo do tempo', 'interpreter', 'latex');
grid minor;

subplot(2,2,4)
plot(t_sim, erro_z2p, 'k--', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro', 'Interpreter', 'latex');
title('Erro do estado $\dot{z_{2}}$ estimado ao longo do tempo', 'interpreter', 'latex');
grid minor;


