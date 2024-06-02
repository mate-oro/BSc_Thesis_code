% folder path for functions
addpath('methods/');

% rhs f and its derivatves
d = -999;
f = @(t, y) d*y^(3);
fd = @(t, y) 3*d*(y^2)*f(t,y);
fy = @(t, y) 3*d*(y^2);

% initial setup
a = 0;
b = 0.5;
N = 10;
y0 = 1;

% solving the IVP
Y1 = LENM2(y0, a, b, N, f, fd, fy, 0.6);
Y2 = theta(f,a,b,y0,N, 1/2);
YA = ones(1,N+1)./sqrt(1 - 2*d*linspace(a,b,N+1));

% calculating the errors
E1 = abs(Y1-YA);
E2 = abs(Y2-YA);
disp("h = " + (b-a)/N + ", LENM2 max: " + max(E1) + ", LENM2 end: " + E1(end) + ", Gauss2 max: " + max(E2) + ", Gauss2 end: " + E2(end))

% plotting the results
figure();
hold on
plot(linspace(a,b,N+1),Y1, "-.");
plot(linspace(a,b,N+1),Y2, "--");
plot(linspace(a,b,100),ones(1,100)./sqrt(1 - 2*d*linspace(a,b,100)),":");
ylim([-0.1 1]);
legend({'LENM2', 'Gauss2', 'Analytic'},'Location','northeast');

% modified iteration for large stepsize
function y = theta(f,a,b,y0,N,r)
    dt = (b-a)/N;
    y = zeros(length(y0),N+1);
    y(:,1) = y0;
    t = a;
    for i = 1:N
        y(:,i+1) = y(:,i) + theta_step(f,y(i),dt,t,r); 
        t = t+dt;
    end
end

function y_step = theta_step(f,y0,dt,t0,r)
    options = optimset('Display','off');
    k1 = f(t0,y0);

    error = @(k) k-f(t0 + r*dt,y0 + r*dt*k);
    
    k1 = fsolve(error, k1, options);

    % Starting point for the calculation of the plot
    %k1_guess = f(t0 + r*dt,y0 + r*dt*k1);
    %k1 = fsolve(error, k1_guess, options);

    y_step = r*dt*k1;
end
