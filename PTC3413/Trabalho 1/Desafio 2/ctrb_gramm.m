function out = ctrb_gramm(A,B,Tf)
f = @(tau) expm(tau*A)*B*B'*expm(tau*A');
Wc = @(T) integral(f, 0, T, 'ArrayValued',1); 
out = Wc(Tf);
end


