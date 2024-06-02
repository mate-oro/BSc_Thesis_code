function y = JE2(f,a,b,y0,N)
    h = (b-a)/N;
    t = linspace(a,b,N+1);
    s = length(y0);
    y = zeros(s,N+1);

    y(:,1) = y0;
    for j=1:N
        k1 = f(t(j),y(:,j));
        y(:,j+1) = y(:,j) + h*f(t(j)+0.5*h,y(:,j)+0.5*h*k1);
    end
end