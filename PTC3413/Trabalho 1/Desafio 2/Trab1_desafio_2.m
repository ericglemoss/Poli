format short;
close all;
clc;
global Tf z0_c zf_c A11 B1
%% Definição dos parâmetros que serão utilizados
M = 4; K = 9; %sqrt(K/M) = 1.5
Tf = 5;
T_u = 0:.01:Tf;
t = 0:.01:10;
%% Condições inciais e finais dos estados x
x10 = 3; x20 = 1;
x10_dot = 0; x20_dot = 0;        

x1f = 3; x2f = 7;
x1f_dot = 1; x2f_dot = 0;
%% Declaração das matrizes originais do problema
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0];

x0 = [x10; x10_dot; x20; x20_dot];
xf = [x1f; x1f_dot; x2f; x2f_dot];

%% Matriz de transformação de base do problema 2: x = Hz => z = H_inv*x
H_inv = [0.5 0 0.5 0; 0 0.5 0 0.5; 0.5 0 -0.5 0; 0 0.5 0 -0.5];
H = inv(H_inv);

%% Mudando a base do problema
A_1 = H_inv*A*H;
B_1 = H_inv*(B*[1;1]);
C_1 = C*H;
z0 = H_inv*x0;
zf = H_inv*xf;

%% Decomposição em parte controlável e não controlável
Co = ctrb(A_1,B_1);
k = rank(Co);
[U,S,V] = svd(Co);
T = U;

A_til = inv(T)*A_1*T;
B_til = inv(T)*B_1;
C_til = C_1*T;

A11 = A_til(1:k,1:k);
A12 = A_til(1:k,k+1:end);
B1 = B_til(1:k,1);
C1 = C_til(:,1:k);


A22 = A_til(k+1:end,k+1:end);
B2 = B_til(k+1:end,k+1:end);
C2 = C_til(:,k+1:end);

%% Cálculo do novo grammiano
W_c = ctrb_gramm(A11,B1,Tf);


%% Definindo os estados iniciais e finais da parte controlável e ñ. controlável
z0_c = z0(1:2,1);
zf_c = zf(1:2,1);

z0_uc = z0(3:4,1);
zf_uc = zf(3:4,1);

%% Declaração dos novos espaços de estados (parte controlávele ñ controlável)
sys_c = ss(A11,B1,C1,[],'StateName',{'z1','z1_dot'});
sys_uc = ss(A22,B2,C2,[],'StateName',{'z2','z2_dot'});

%% Cálculo da entrada necessária para mudança de estado
 u = zeros(length(t),1);

for k = 1:length(T_u)
  u(k,1) = -sys_c.B'*expm((Tf-T_u(k))*sys_c.A')*inv(W_c)*(expm(Tf*sys_c.A)*z0_c-zf_c);
end

%% Simulação da parte controlável via lsim
[yc,tc,zc] = lsim(sys_c,u,t,z0_c);

figure
subplot(2,2,1)
plot(tc, zc(:,1), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l = legend('$z_{1}$');
set(l,'Interpreter', 'latex');
grid minor;

subplot(2,2,2)
plot(tc, zc(:,2), 'r--', 'LineWidth', 1);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l = legend('$\dot{z_{1}}$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,2,[3,4])
plot(tc, yc);
xlabel('Tempo (s)');
ylabel('Valor do estado');
l = legend('$y_{1}$', '$y_{2}$');
set(l, 'Interpreter', 'latex');
grid minor;


%% Simulação do modelo via Simulink e disposição das saídas do mesmo
sim('Trab1_desafio2',10);

%% Demonstração dos resultados
figure('Name','Estados simulados pelo Simulink')
subplot(2,2,1)
plot(t_sim, z1_sim(:,1), 'r', t_sim, z1_sim(:,2), 'k','LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Vetor de estados $z_{1}^{1}$ ao longo do tempo', 'Interpreter', 'latex');
l = legend('$z_{1}$','$\dot{z_{1}}$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,2,2)
plot(t_sim, z1ponto_sim(:,1), 'r', t_sim, z1ponto_sim(:,2), 'k','LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Vetor de estados $\dot{z_{1}}^{1}$ ao longo do tempo', 'Interpreter', 'latex');
l = legend('$\dot{z_{1}}$','$\ddot{z_{1}}$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,2,3)
plot(t_sim, z2_sim(:,1), 'c', t_sim, z2_sim(:,2), 'b','LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Vetor de estados $z_{2}^{2}$ ao longo do tempo', 'Interpreter', 'latex');
l = legend('$z_{2}$','$\dot{z_{2}}$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,2,4)
plot(t_sim, z2ponto_sim(:,1), 'c', t_sim, z2ponto_sim(:,2), 'b', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor do estado', 'Interpreter', 'latex');
title('Vetor de estados $\dot{z_{2}}^{2}$ ao longo do tempo', 'Interpreter', 'latex');
l = legend('$\dot{z_{2}}$','$\ddot{z_{2}}$');
set(l, 'Interpreter', 'latex');
grid minor;

figure('Name', 'Entrada simulada via Simulink');
plot(t_sim, u_sim, 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor da entrada', 'Interpreter', 'latex');
title('Entrada calculada ao longo do tempo', 'Interpreter', 'latex');
l = legend('$u$');
set(l, 'Interpreter', 'latex');
grid minor;

figure('Name', 'Saídas simuladas via simulink');
subplot(2,1,1);
plot(t_sim, y1_sim(:,1), 'r', t_sim, y1_sim(:,2), 'k', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor da saida', 'Interpreter', 'latex');
title('Vetor de saidas $y_{1}^{1}$ da parte controlavel', 'interpreter', 'latex');
l = legend('$y_{1}$', '$y_{2}$');
set(l, 'Interpreter', 'latex');
grid minor;

subplot(2,1,2)
plot(t_sim, y2_sim(:,1), 'c', t_sim, y2_sim(:,2), 'm', 'LineWidth', 1);
xlabel('Tempo (s)', 'Interpreter', 'latex');
ylabel('Valor da saida', 'Interpreter', 'latex');
title('Vetor de saidas $y_{2}^{2}$ da parte nao controlavel', 'interpreter', 'latex');
l = legend('$y_{1}$','$y_{2}$');
set(l, 'Interpreter', 'latex');
grid minor;





