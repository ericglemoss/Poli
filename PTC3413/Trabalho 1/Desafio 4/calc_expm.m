function out = calc_expm(u)
global A22 C2
out = expm(u(2)*A22')*C2'*u(1);
end

