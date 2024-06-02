function [Steps, y]=DOPRI5(f,t0,T,y0,epsnull)
    s = length(y0);
    h = sqrt(epsnull);
    minstep = 1e-20; 
    maxit=1e8;
    y=zeros(length(y0),maxit);
    Steps = zeros(1,maxit);
    Steps(1) = t0;
    a      =  zeros(7,7);
    a(2,1) =  0.2;
    a(3,1) =  3/40;
    a(3,2) =  9/40;
    a(4,1) =  44/45;
    a(4,2) = -56/15;
    a(4,3) =  32/9;
    a(5,1) =  19372/6561;
    a(5,2) = -25360/2187;
    a(5,3) =  64448/6561;
    a(5,4) = -212/729;
    a(6,1) =  9017/3168;
    a(6,2) = -355/33;
    a(6,3) =  46732/5247;
    a(6,4) =  49/176;
    a(6,5) = -5103/18656;
    a(7,1) =  35/384;
    a(7,3) =  500/1113;
    a(7,4) =  125/192;
    a(7,5) = -2187/6784;
    a(7,6) =  11/84;

    c     = zeros(1,7); 
    c(1)  = 0;
    c(2)  = 0.2;
    c(3)  = 0.3; 
    c(4)  = 0.8;
    c(5)  = 8/9;
    c(6)  = 1;
    c(7)  = 1;
    
    er = zeros(1,7);
    er(1) =  71/57600;
    er(2) =  0.0;
    er(3) = -71/16695;
    er(4) =  71/1920;
    er(5) = -17253/339200;
    er(6) =  22/525;
    er(7) = -1/40;

    b = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84; 0];

    y(:,1)=y0;
    t = t0;
    j=1;
    hmin = h;
    k = zeros(s,7);
    while t<T && j<maxit
        for i=1:7 
            k(:,i) = f(t + c(i)*h, y(:,j) + h*k*(a(i,:)'));
        end
        yerr = h*norm(k*er');
        if yerr <= epsnull
            y(:,j+1) = y(:,j) + h*k*b;
            j = j+1;
            t = t+h;
            Steps(j) = t;
        elseif h == minstep
            disp("Too samll step required at t = " + t );
            break
        end
        if yerr ~=0
            fac = max(min((epsnull/yerr)^(1/5),2),1e-1);
            h = 0.9*h*fac;
            h = max(h,minstep);
            hmin = min(h, hmin);
            h = min(h,T-t);
        else
            h = 2*h;
        end
    end
    y = y(:,1:j);
    Steps = Steps(1:j);
    
    disp("Minimum h = " + num2str(hmin))
    
    if j == maxit
        error('Max iteration reached at t = ' + num2str(t));
    end
end