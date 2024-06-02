function y_step = LobIIID_step(f,y0,dt,t0)
    options = optimset('Display','off');
    s = length(y0);

    a11 = 1/2;
    a12 = 1/2;
    a21 = -1/2;
    a22 = 1/2;

    b1 = 1/2;
    b2 = 1/2;

    c1 = 0;
    c2 = 1;

    k=f(t0,y0);

    error=@(k) [k(1:s) - f(t0 + c1*dt,y0 + dt*(a11*k(1:s)+a12*k(s+1:2*s)));
                k(s+1:2*s) - f(t0 + c2*dt,y0 + dt*(a21*k(1:s) + a22*k(s+1:2*s)))];
    
    [k] = fsolve(error, [k;k], options);
    y_step = dt*(b1*k(1:s) + b2*k(s+1:2*s));
end