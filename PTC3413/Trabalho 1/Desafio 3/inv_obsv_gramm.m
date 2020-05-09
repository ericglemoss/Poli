function out = inv_obsv_gramm()
global A C T
g = @(tau) expm(tau*A')*C'*C*expm(tau*A);
Wobs = @(Tf) integral(g, 0, Tf, 'ArrayValued',1);
out = inv(Wobs(T));
end

