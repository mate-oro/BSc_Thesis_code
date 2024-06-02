function y = RK4(f,a,b,y0,N)
    h = (b-a)/N;
    t = linspace(a,b,N+1);
    s = length(y0);
    y = zeros(s,N+1);

    y(:,1) = y0;
    for j = 1:N
        k1 = f(t(j),y(:,j));
        k2 = f(t(j)+0.5*h,y(:,j)+0.5*h*k1);
        k3 = f(t(j)+0.5*h,y(:,j)+0.5*h*k2);
        k4 = f(t(j)+h,y(:,j)+h*k3);
        y(:,j+1) = y(:,j)+h*(1/6*k1+1/3*k2+1/3*k3+1/6*k4);
    end
end