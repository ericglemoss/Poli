function out = obsv_gramm(T)
%     f = @(tau) expm(tau*A.')*C.'*C*expm(tau*A);           
%     W = @(t) integral(f, 0, t, 'ArrayValued',1);        
%     out = W(T)
      opt = gramOptions('TimeIntervals', [0 T]);
      out = gram(sys,'o',opt);
end

