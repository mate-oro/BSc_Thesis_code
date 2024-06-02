function y = AENM2(y0, a, b, N, f, fd)
    dt = (b-a)/N;
    t = linspace(a,b,N+1);
    y = zeros(length(y0),N+1);
    y(:,1) = y0;
    for i=1:N
        y(:,i+1) = y(:,i) + AENM2_step(y(:,i), t(i), dt, f, fd);
    end
end