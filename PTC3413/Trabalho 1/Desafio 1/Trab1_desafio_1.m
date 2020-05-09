format short;
close all;
clc;
global A B C Tf x0 xf
%% Definição dos parâmetros que serão utilizados
t = 0:.1:10; % Tempo de simulação 
M = 4; K = 9; %sqrt(K/M) = 1.5
Tf = 5; % Tempo no qual a entrada age
T = 0:.1:Tf;

%% Condições iniciais e finais dos estados de x
x10 = 1; x20 = 1;
x1f = 6; x2f = 6;

%% Definição das matrizes do espaço de estados e das condições iniciais e finais
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0];

x0 = [x10; 0; x20; 0];
xf = [x1f; 0; x2f; 0];

sys = ss(A,B,C,[],'StateName',{'x1','x1_dot','x2','x2_dot'});

%% Checando se o sistema é controlável
if rank(ctrb(sys.A,sys.B)) == 4
    disp('Sistema do Desafio 1 é controlável');
    Co = ctrb(sys.A,sys.B);
else
    error('Sistema do Desafio 1 não é controlável');
end

%% Cálculo das entradas necessárias para que em um tempo Tf x1(Tf) = x2(Tf) = xf
Wc = ctrb_gramm(); % Grammiano de controlabilidade

u_1 = zeros(length(t),1);
u_2 = zeros(length(t),1);

for k = 1:length(T)
    tmp = -sys.B'*expm((Tf-T(k))*sys.A')*inv(Wc)*(expm(Tf*sys.A)*x0-xf);
    u_1(k,1) = tmp(1,1);
    u_2(k,1) = tmp(2,1);
end

uf = [u_1 u_2];

[y,t,x] = lsim(sys,uf,t,x0);


figure('Name','Estados obtidos através da simulação via lsim');
subplot(2,2,1)
plot(t, x(:,1), 'r--', t, x(:,2), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l_1 = legend('$x_{1}$','$\dot{x_{1}}$');
set(l_1,'Interpreter','latex');
grid minor;

subplot(2,2,2)
plot(t, x(:,3), 'r--', t, x(:,4), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l_2 = legend('$x_{2}$','$\dot{x_{2}}$');
set(l_2,'Interpreter','latex');
grid minor;

subplot(2,2,[3,4])
plot(t, y(:,1), 'b--', t, y(:,2), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l_3 = legend('$y_{1}$', '$y_{2}$');
set(l_3, 'Interpreter', 'latex');
grid minor;

%% Verificando se o problema é controlável apenas para a entrada u1
if rank(ctrb(sys.A,sys.B(:,1))) == 4
    disp('O sistema do Desafio 1 é controlável apenas para u1.');
else
    disp('O sistema do Desafio 1 não é controlável apenas para u1.');
end

%% Simulação do modelo no Simulink
 sim('Trab1_desafio1',10);
 
%% 
figure('Name','Estados obtidos a partir da simulação via Simulink');
subplot(2,2,1);
plot(t_sim, x1_t, 'r', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
title('Valor do estado $x_1$ ao longo do tempo', 'Interpreter','latex');
l = legend('$x_{1}$');
set(l,'Interpreter','latex');
grid minor;


subplot(2,2,2)
plot(t_sim, x1p_t, 'k', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
title('Valor do estado $\dot{x_{1}}$ ao longo do tempo', 'Interpreter','latex');
l = legend('$\dot{x_{1}}$');
set(l,'Interpreter','latex');
grid minor;


subplot(2,2,3)
plot(t_sim, x2_t, 'b', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
title('Valor do estado $x_2$ ao longo do tempo', 'Interpreter','latex');
l = legend('$x_{2}$');
set(l,'Interpreter','latex');
grid minor;

subplot(2,2,4)
plot(t_sim, x2p_t, 'm', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
title('Valor do estado $\dot{x_2}$ ao longo do tempo', 'Interpreter','latex');
l = legend('$\dot{x_{2}}$');
set(l,'Interpreter','latex');
grid minor;


figure('Name', 'Saídas calculadas pelo modelo do Simulink');
subplot(2,1,1);
plot(t_sim, y1_t, 'r', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor da saída');
title('Valor do estado $y_1$ ao longo do tempo', 'Interpreter','latex');
l = legend('$y_1$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,1,2);
plot(t_sim, y2_t, 'k', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor da saída');
title('Valor do estado $y_2$ ao longo do tempo', 'Interpreter','latex');
l = legend('$y_2$');
set(l, 'Interpreter', 'latex');
grid minor;

