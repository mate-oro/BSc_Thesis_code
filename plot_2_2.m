% folder path for functions
addpath('methods/');

% rhs f
f = @(t,y) y;

% initial setup
a = 0;
b = 2;
h = 0.1;
%h = 0.25;
%h = 0.5;
N = ceil((b-a)/h);
y0 = 1;

% solving the IVP
EE = Theta(f,a,b,y0,N,0);
JE = JE2(f,a,b,y0,N);

% plotting the results
X = linspace(a,b,N+1);
XA = linspace(a,b,100);
plot(XA,exp(XA),'-')
hold on
plot(X,EE,'*')
plot(X,JE,'o')
ylim([1,8])
legend({'Pontos', 'EE','JE'},'Location','northwest','FontSize',10);
hold off;