function u = calc_u(input)
global A B Tf x0 xf

Wc = ctrb_gramm();

if input < Tf
    u = -B'*expm((Tf-input)*A')*inv(Wc)*(expm(Tf*A)*x0-xf);
else
    u = [0;0];
end


