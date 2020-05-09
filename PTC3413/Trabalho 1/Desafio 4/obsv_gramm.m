function out = obsv_gramm(A,C,Tf)
g = @(tau) expm(tau*A')*C'*C*expm(tau*A);
Wobs = @(Tf) integral(g, 0, Tf, 'ArrayValued',1);
out = Wobs(Tf);
end

