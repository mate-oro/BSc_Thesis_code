function [Steps, y] = AGauss2(f,t0,T,y0,epsnull)
    h = 1e-3;
    r = 1/2;
    minstep = (T-t0)*(1e-20); 
    maxit = 1e8;
    y = zeros(length(y0),maxit);
    Steps = zeros(1,maxit);
    Steps(1) = t0;
    y(:,1) = y0;
    t = t0;
    j = 1;
    hmin = h;
    while t<T && j<maxit
        y1 = y(:,j) + Theta_step(f, y(:,j), h, t, r);
        y21 = y(:,j) + Theta_step(f, y(:,j), h/2, t, r);
        y22 = y21 + Theta_step(f, y21, h/2, t+h/2, r);
        yerr = norm(y1-y22,'inf');
        if yerr <= epsnull
            y(:,j+1) = y22;
            j = j+1;
            t = t+h;
            Steps(j) = t;
        elseif h == minstep
            disp("Too samll step required at t = " + t );
            break
        end
        if yerr ~= 0
            fac = max(min((epsnull/yerr)^(1/2),2),1e-1);
            h = 0.9*h*fac;
            h = max(h,minstep);
            hmin = min(h, hmin);
            h = min(h,T-t);
        end
    end
    y = y(:,1:j);
    Steps = Steps(1:j);
    disp("Minimum h = " + num2str(hmin));
    if j == maxit
        error('Max iteration reached at t = ' + num2str(t));
    end
end