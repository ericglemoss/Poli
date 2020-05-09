function u = calcula_entrada(inp)
global A11 B1 Tf z0_c zf_c
Wc = ctrb_gramm(A11,B1,Tf);
 if inp < Tf 
     u = -B1'*expm((Tf-inp)*A11')*inv(Wc)*(expm(Tf*A11)*z0_c-zf_c);
 else
    u = 0;
end
end
