function xo = calc_xo(u)
    opt = gramOptions('TimeIntervals', [0 u(2)]);
    Wo = gram(sys,'o',opt);
    @w = int

end
