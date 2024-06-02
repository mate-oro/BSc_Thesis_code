function y_step = Theta_step(f,y0,dt,t0,r)
    options = optimset('Display','off');
    k1 = f(t0,y0);
    k1_guess = f(t0 + r*dt,y0+dt*k1);
    
    error = @(k) k-f(t0 + r*dt,y0 + dt*r*k);
    
    k1 = fsolve(error, k1_guess, options);
    y_step = dt*k1;
end