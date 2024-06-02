% folder path for functions
addpath('methods/');

% rhs f
f = @(t,y) -9*y;

% initial setup
a = 0;
b = 2;
h = 0.1;
%h = 0.25;
%h = 0.5;
N = ceil((b-a)/h);
y0 = 1;

% solving the IVP
IE = Theta(f,a,b,y0,N,1);
JE = JE2(f,a,b,y0,N);
T = Theta(f,a,b,y0,N,1/4);
RK = RK4(f,a,b,y0,N);

% plotting the results
X = linspace(a,b,N+1);
XA = linspace(a,b,100);
plot(XA,exp(-9*XA),'--')
hold on
plot(X,JE,'-')
plot(X,T,'.-')
plot(X,RK,'*')
plot(X,IE,'o')
ylim([-0.5,2])
legend({'Pontos','JE','Theta 1/4','RK4','IE'},'Location','northeast')
hold off;