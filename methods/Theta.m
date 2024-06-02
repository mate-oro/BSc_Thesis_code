function y = Theta(f,a,b,y0,N,r)
    dt = (b-a)/N;
    y = zeros(length(y0),N+1);
    y(:,1) = y0;
    t = a;
    for i = 1:N
        y(:,i+1) = y(:,i) + Theta_step(f,y(:,i),dt,t,r); 
        t = t+dt;
    end
end
