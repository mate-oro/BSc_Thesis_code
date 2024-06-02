function y = AENM2_step(y0, t, dt, f, fd)
    f0 = f(t,y0);
    fd0 = fd(t,y0);
    y = (2*dt*(f0.^2))./(2*f0 - dt*fd0);
end