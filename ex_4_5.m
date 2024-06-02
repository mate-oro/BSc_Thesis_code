% folder path for functions
addpath('methods/');

% rhs f and its derivatves
f = @(t, y) y^2 - exp(-2000*t) - 1002*exp(-1000*t)-1;
fd = @(t, y) 2000*exp(-2000*t) + 1002000*exp(-1000*t) + (2*y)*f(t,y);
fy = @(t, y) 2*0;

% initial setup
a = 0;
b = 0.1;
N = 100;
y0 = 2;

% solving the IVP
Y1 = LENM2(y0, a, b, N, f, fd, fy, 0.55);
Y2 = AENM2(y0, a, b, N, f, fd);

% calculating the errors
X = linspace(a,b,N+1);
E1 = abs(Y1-exp(-1000*X)-1);
E2 = abs(Y2-exp(-1000*X)-1);
disp("h = " + (b-a)/N + ", LENM2 max: " + max(E1) + ", LENM2 end: " + E1(end) + ", AENM2 max: " + max(E2) + ", AENM2 end: " + E2(end))

% plotting the results
figure();
hold on;
plot(X,E1, "-.");
plot(X,E2, "-.");
legend({'LENM2', 'AENM2'},'Location','northeast');
hold off;