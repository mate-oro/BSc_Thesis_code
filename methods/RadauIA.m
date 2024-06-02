function y = RadauIA(f,a,b,y0,N)
    h = (b-a)/N;
    y = zeros(length(y0),N+1);
    t = linspace(a,b,N+1);
    y(:,1) = y0;
    for i=1:N
        y(:,i+1) = y(:,i) + RadauIA_step(f,y(:,i),h,t(i));
    end
end