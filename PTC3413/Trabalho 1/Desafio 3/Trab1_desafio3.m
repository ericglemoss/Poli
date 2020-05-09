clc;
%% Inicializações
global A C T
A = [0 1 0 0; -2.25 0 2.25 0; 0 0 0 1; 2.25 0 -2.25 0];
C = [1 0 0 0; 0 0 1 0];
T = 8;
inv_Wo = inv_obsv_gramm();
x0 = [0;1;0;-1];
%% Simulação
sim('Trab1_desafio_3',10);

%% Plot dos estados estimados
figure('Name','Estados iniciais determinados via Simulink');
subplot(2,2,1)
plot(t_sim, x1_calc, 'r', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Estado $x_{1}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;


subplot(2,2,2)
plot(t_sim, x1ponto_calc, 'k', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Estado $\dot{x_{1}}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;

subplot(2,2,3)
plot(t_sim, x2_calc, 'c', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Estado $x_{2}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;

subplot(2,2,4)
plot(t_sim, x2ponto_calc, 'm', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Estado $\dot{x_{2}}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;

%% Plot dos erros
figure('Name','Erros relativos na estimação de cada estado');
subplot(2,2,1)
plot(t_sim,x0(1,1)-x1_calc, 'r', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro relativo', 'Interpreter', 'latex');
title('Erro relativo de $x_{1}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;

subplot(2,2,2)
plot(t_sim,x0(2,1)-x1ponto_calc, 'k', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro relativo', 'Interpreter', 'latex');
title('Erro relativo de $\dot{x_{1}}$ estimado ao longo do tempo' ,'Interpreter', 'latex');
grid minor;

subplot(2,2,3)
plot(t_sim,x0(3,1)-x2_calc, 'c', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro relativo', 'Interpreter', 'latex');
title('Erro relativo de $x_{2}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;

subplot(2,2,4)
plot(t_sim,x0(4,1)-x2ponto_calc, 'm', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do erro relativo', 'Interpreter', 'latex');
title('Erro relativo de $\dot{x_{2}}$ estimado ao longo do tempo', 'Interpreter', 'latex');
grid minor;
