function u = calc_u(input)
A = [0 1 0 0; -1.5 0 1.5 0; 0 0 0 1; 1.5 0 -1.5 0];
B = [0 0; 1.5 0; 0 0; 0 1.5];

x0 = [1; 0; 1; 0];
xf = [7; 0; 7; 0];

Wc = ctrb_gramm(A,B,input(2))
u = -B'*expm((input(2)-input(1))*A')*inv(Wc)*(expm(input(2)*A)*x0-xf);

end

