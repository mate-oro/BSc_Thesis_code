function y = LENM2_step(y0, t0, dt, f, fd, fy, a)
    yd = diag(y0);
    f0 = f(t0,y0);
    fd0 = fd(t0,y0);
    fy0 = fy(t0,y0);
    y = 2*(yd*y0 + dt*yd*f0 - dt*a*yd^2*fy0)./(2*y0 - 2*dt*a*yd*fy0 - dt^2*fd0 + 2*dt^2*a*fy0.*f0);
end