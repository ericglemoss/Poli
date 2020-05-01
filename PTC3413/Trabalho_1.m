%% Definição das variáveis que serão utilizadas
% syms t
M = 4; K = 9; Tf = 5; x1f = 15; x2f = 15; t = 15;
format short; 
close all;

%% Definição das matrizes do desafio 1 e do respec. espaço de estados
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0]; 
x0 = [10; 0; 10; 0];
sys = ss(A,B,C,[], 'StateName',{'x1','x1_d','x2','x2_d'});

%% Checa se o sistema é controlável
if rank(ctrb(sys.A,sys.B)) == 4
    disp('Sistema do desafio 1 é controlável')
    Co = ctrb(A,B);
else
    error('Sistema do desafio 1 não é controlável')
end

%% Cálculo das entradas necessárias para que em um tempo T x1(T) = x2(T) = xf
xf = [x1f; 0; x2f; 0];
Wc = ctrb_gramm(A,B,Tf);
uf = -sys.B'*expm((Tf-t)*A')*inv(Wc)*(expm(Tf*A)*x0-x1f);

%% Verificando se o problema é resolvível apenas para a entrada u1
if rank(ctrb(sys.A,sys.B(:,1))) == 4
    disp('O sistema do desafio 1 é controlável apenas para u1.');
else
    disp('O sistema do desafio 1 não é controlável apenas para u1.');
end

%% Desafio 2. Reescrevendo o modelo (A_1,B_1) e declaração das novas variáveis(Item A)
T_inv = [0.5 0 0.5 0; 0 0.5 0 0.5; 0.5 0 -0.5 0; 0 0.5 0 -0.5]; % x = Tz => z = T_inv*x;
T = inv(T_inv);
z0 = [10;0;10;0];
zf = [15;0;15;0];

A_1 = T_inv*A*T;
tmp = B(:,1)+B(:,2); B_1 = T_inv*tmp;
C_1 = C*T;

%% Mostrando que a dinâmica (z1,z1_dot) corresponde a parte controlável e que (z2,z2_dot) a não controlável
[Abar,Bbar,Cbar] = ctrbf(A_1,B_1,C_1);

t_1 = 0:.01:15;

Co_1 = ctrb(A_1,B_1);
rank_co1 = rank(Co_1);

A_c = Abar(rank_co1+1:end,rank_co1+1:end);
A_uc = Abar(1:rank_co1,1:rank_co1);

%% Cálculo dos grammianos do desafio 2
Wc_bar = ctrb_gramm(Abar,Bbar,Tf);
Wo_bar = obsv_gramm(Abar,Cbar,Tf);

%% Declaração do novo espaço de estados
sys_bar = ss(Abar,Bbar,Cbar,[],'StateName',{'z1','z1_dot','z2','z2_dot'});

%% Cálculo da entrada para levar o sistema de z1(0) até z1(Tf)
u_zf = sys_bar.B'*expm((Tf-t)*sys_bar.A')*inv(Wo_bar)*(expm(Tf*sys_bar.A)*z0-zf);
teste = zeros(length(t_1),1);
teste(1,1) = u_zf;

[y_bar,t_bar,x_bar] = lsim(sys_bar,teste,t_1,z0);
plot(t_bar,x_bar);












 
    


