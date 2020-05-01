%% Definição das variáveis que serão utilizadas
% syms t
M = 4; K = 9; Tf = 5; x1f = 10; x2f = 10; t = 15;
format short; 

%% Definição das matrizes do desafio 1 e do respec. espaço de estados
A = [0 1 0 0; -K/M 0 K/M 0; 0 0 0 1; K/M 0 -K/M 0];
B = [0 0; 1/M 0; 0 0; 0 1/M];
C = [1 0 0 0; 0 0 1 0]; 
x0 = [10; 0; 10; 0];
sys = ss(A,B,C,[], 'StateName',{'x1','x1_d','x2','x2_d'});
isstable(sys);

%% Opções para o cálculo dos Grammianos
opt = gramOptions('TimeIntervals', [0 Tf]);

%% Checa se o sistema é controlável
if rank(ctrb(sys.A,sys.B)) == 4
    disp('Sistema do desafio 1 é controlável')
    Co = ctrb(A,B);
else
    error('Sistema do desafio 1 não é controlável')
end

%% Cálculo dos Grammianos de Controlabilidade e Observabilidade
Wc = gram(sys,'c',opt);
Wo = gram(sys,'o',opt);

%% Cálculo das entradas necessárias para que em um tempo T x1(T) = x2(T) = xf
xf = [x1f; 0; x2f; 0];
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

A_1 = T_inv*A*T;
tmp = B(:,1)+B(:,2); B_1 = T_inv*tmp;
C_1 = C*T;
z = [z1; z1_dot; z2; z2_dot];

%% Mostrando que a dinâmica (z1,z1_dot) corresponde a parte controlável e que (z2,z2_dot) a não controlável
[Abar,Bbar,Cbar] = ctrbf(A_1,B_1,C_1);
Co_1 = ctrb(A_1,B_1);
rank_co1 = rank(Co_1);
A_c = Abar(rank_co1+1:end,rank_co1+1:end);
A_uc = Abar(1:rank_co1,1:rank_co1);
