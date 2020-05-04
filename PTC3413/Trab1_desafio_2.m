format short;
close all;
t = 0:.01:15;
%% Definição dos parâmetros que serão utilizados
M = 4; K = 9; %sqrt(K/M) = 1.5
Tf = 5;

%% Condições inciais e finais dos estados x
x10 = 2; x20 = 2;
x1f = 10; x2f = 10;

%% Declaração das matrizes originais do problema
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0];

x0 = [x10; 0; x20; 0];
xf = [x1f; 0; x2f; 0];

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

A_c = A_til(1:k,1:k);
B_c = B_til(1:k,1);
C_c = C_til(:,1:k);

A_uc = A_til(k+1:end,k+1:end);
B_uc = B_til(k+1:end,k+1:end);
C_uc = C_til(:,k+1:end);

%% Cálculo dos novos grammianos
W_c = ctrb_gramm(A_c,B_c,Tf);
W_uc = ctrb_gramm(A_uc,B_uc,Tf);

z0_c = z0(1:2,1);
zf_c = zf(1:2,1);

z0_uc = z0(3:4,1);
zf_uc = zf(3:4,1);

%% Declaração dos novos espaços de estados (parte controlávele ñ controlável)
sys_c = ss(A_c,B_c,C_c,[],'StateName',{'z1','z1_dot'});
sys_uc = ss(A_uc,B_uc,C_uc,[],'StateName',{'z2','z2_dot'});

%% Cálculo da entrada necessária para mudança de estado
for k = 1:length(t)
  u(k,1) = -sys_c.B'*expm((Tf-t(k))*sys_c.A')*inv(W_c)*(expm(Tf*sys_c.A)*z0_c-zf_c);
end

%% Simulação da parte controlável
[yc,tc,xc] = lsim(sys_c,u,t,z0_c);
plot(tc,xc);



