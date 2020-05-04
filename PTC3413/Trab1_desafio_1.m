format short;
close all;
t = 0:.1:10;
%% Definição dos parâmetros que serão utilizados
M = 4; K = 9; %sqrt(K/M) = 1.5
Tf = 5;

%% Condições iniciais e finais dos estados de x
x10 = 1; x20 = 1;
x1f = 10; x2f = 10;

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
Wc = ctrb_gramm(sys.A,sys.B,Tf); % Grammiano de controlabilidade

u_1 = zeros(length(t),1);
u_2 = zeros(length(t),1);

for k = 1:length(t)
    tmp = -sys.B'*expm((Tf-t(k))*sys.A')*inv(Wc)*(expm(Tf*sys.A)*x0-xf);
    u_1(k,1) = tmp(1,1);
    u_2(k,1) = tmp(2,1);
end

uf = [u_1 u_2];

[y,t,x] = lsim(sys,uf,t,x0);
plot(t, x(:,1), 'r--', t, x(:,2), 'm--', 'Linewidth', 1);
l = legend('$x_{1}$','$x_{2}$');
set(l,'Interpreter','latex');

%% Verificando se o problema é controlável apenas para a entrada u1
if rank(ctrb(sys.A,sys.B(:,1))) == 4
    disp('O sistema do Desafio 1 é controlável apenas para u1.');
else
    disp('O sistema do Desafio 1 não é controlável apenas para u1.');
end
 
% sim('Trab1_desafio1',10);