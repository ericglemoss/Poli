function out = calc_expm(u)
global A C
out = expm(u(3,1)*A')*C'*u(1:2,1);
end

