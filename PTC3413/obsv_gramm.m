function out = obsv_gramm(A,C,T)
g = @(tau) expm(tau*A')*C'*C*expm(tau*A);
Wobs = @(T) integral(g, 0, T, 'ArrayValued',1);
out = Wobs(T);
end

