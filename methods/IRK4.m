function y = IRK4(f,a,b,y0,N)
    h = (b-a)/N;
    y = zeros(length(y0),N+1);
    y(:,1) = y0;
    t = a;
    for i=1:N
        y(:,i+1) = y(:,i) + IRK4_step(f,y(:,i),h,t);
        t = t+h;
    end
end