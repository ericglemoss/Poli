function out = obsv_gramm()
global A C Tf
g = @(tau) expm(tau*A')*C'*C*expm(tau*A);
Wobs = @(T) integral(g, 0, T, 'ArrayValued',1);
out = Wobs(Tf);
end

