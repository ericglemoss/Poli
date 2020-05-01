function out = ctrb_gramm(A,B,T)
f = @(tau) expm(A*tau)*B*B'*expm(A'*tau);
W = @(T) integral(f, 0, T, 'ArrayValued',1); 
out = W(T);
end


