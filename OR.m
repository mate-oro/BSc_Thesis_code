% folder path for functions
addpath('methods/')

% initial setup
a = 0;
b = 1;
y0 = [1; 1];
N1 = 10;
N2 = 10*N1;
X1 = linspace(a,b,N1+1);
X2 = linspace(a,b,N2+1);

% setup of the values for which we calculete the order 
M = 50;
E = 5.9;
lam = linspace(-1,E,M);
lam = 10.^lam;
e1 = zeros(1,M);
e2 = zeros(1,M);
e3 = zeros(1,M);
e4 = zeros(1,M);
e5 = zeros(1,M);
E1 = exp([-2*X1; -X1]);
E2 = exp([-2*X2; -X2]);
for i=1:M
    f = @(t,y) ([-(lam(i) + 2)*y(1) + lam(i)*y(2)^2; y(1) - y(2) - y(2)^2]);
    T1 = IRK4(f,a,b,y0,N1);
    T2 = IRK4(f,a,b,y0,N2);
    e1(i) = log10(max(norm(T1-E1, 1),[],'all')/max(norm(T2-E2, 1),[],'all'));
    T1 = theta(f,a,b,y0,N1,1);
    T2 = theta(f,a,b,y0,N2,1);
    e5(i) = log10(max(norm(T1-E1, 1),[],'all')/max(norm(T2-E2, 1),[],'all'));
    T1 = RadauIIA(f,a,b,y0,N1);
    T2 = RadauIIA(f,a,b,y0,N2);
    e2(i) = log10(max(norm(T1-E1, 1),[],'all')/max(norm(T2-E2, 1),[],'all'));
    T1 = LobIIID(f,a,b,y0,N1);
    T2 = LobIIID(f,a,b,y0,N2);
    e3(i) = log10(max(norm(T1-E1, 1),[],'all')/max(norm(T2-E2, 1),[],'all'));
    disp(i); % Just to indicate how far the calculation has progressed
end

% for some methods the standard grid is not fine enough to converge in the very stiff
% case, so wee choose some factor MM times finer computational grid 
MM = 2.3;
N10 = MM*N1;
N20 = 10*N10;
X10 = linspace(a,b,N10+1);
X20 = linspace(a,b,N20+1);
E10 = exp([-2*X10; -X10]);
E20 = exp([-2*X20; -X20]);
lam10 = linspace(-1,E,MM*M);
lam10 = 10.^lam10;
for i=1:2*M
   f = @(t,y) ([-(lam10(i) + 2)*y(1) + lam10(i)*y(2)^2; y(1) - y(2) - y(2)^2]);
   T1 = RadauIA(f,a,b,y0,N10);
   T2 = RadauIA(f,a,b,y0,N20);
   e4(i) = log10(max(norm(T1-E10, 1),[],'all')/max(norm(T2-E20, 1),[],'all'));
   disp(i); % Just to indicate how far the calculation has progressed
end
for i=2*M+1:MM*M
    e4(i) = e4(i-1);
end

% plotting the results
figure()
subplot(3,1,[1 2]);
hold on
plot(lam,e1,'-', 'Color','#00cd6c','LineWidth', 2,'MarkerSize',5);
plot(lam,e2,'-', 'Color','#af58ba','LineWidth', 2,'MarkerSize',5);
plot(lam10,e4,'-', 'Color','#009ade','LineWidth', 2,'MarkerSize',5);
plot(lam,e3,'-', 'Color','#ffc61e','LineWidth', 2,'MarkerSize',5);
plot(lam,e5,'-', 'Color','#f28522','LineWidth', 2,'MarkerSize',5);
set(gca,'XScale','log');
legend({'Gauss4','RadauIIA','RadauIA','LobattoIIID','IE'},'Location','northeast','FontSize',12)
hold off

% modified iteration for the more stiff case
function y = theta(f,a,b,y0,N,r)
    dt = (b-a)/N;
    y = zeros(length(y0),N+1);
    y(:,1) = y0;
    t = a;
    for i = 1:N
        y(:,i+1) = y(:,i) + theta_step(f,y(:,i),dt,t,r); 
        t = t+dt;
    end
end

function y_step = theta_step(f,y0,dt,t0,r)
    options = optimset('Display','off');
    k1 = f(t0,y0);
    
    error = @(k) k-f(t0 + r*dt,y0 + dt*r*k);

    % Here we dont use the gess as it seems to slow down convergence
    %k1_guess = f(t0 + r*dt,y0 + r*dt*k1);
    %k1 = fsolve(error, k1_guess, options);

    k1 = fsolve(error, k1, options);
    y_step = r*dt*k1;
end
