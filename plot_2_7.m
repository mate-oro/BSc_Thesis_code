% folder path for functions
addpath('methods/');

% rhs f and its derivatves
f = @(t,y) (-999*(y-cos(t)));

% initial setup
a = 0;
b = 1.5;
y0 = 0;
N = 30;

% solving the IVP
T1 = Theta(f,a,b,y0,60,1/2);
T2 = Theta(f,a,b,y0,30,1/2);
IE = Theta(f,a,b,y0,N,1);

% plotting the results
X = linspace(a,b,N+1);
XA = linspace(a,b,10*N+1);
figure()
subplot(3,1,[1 2]);
hold on
plot(linspace(a,b,30+1),T2(1,:),'o-', 'Color','#00cd6c','MarkerFaceColor','#00cd6c','LineWidth', 1.25,'MarkerSize',5)
plot(linspace(a,b,60+1),T1(1,:),'o-', 'Color','#af58ba','MarkerFaceColor','#af58ba','LineWidth', 1.25,'MarkerSize',5)
plot(X,IE(1,:),'o-','Color','#ffc61e','LineWidth', 1.5,'MarkerFaceColor','white','MarkerSize',7)
ylim([0,1.95])
legend({'Gauss2 30 lépés','Gauss2 60 lépés','IE 30 lépés'},'Location','northeast','FontSize',12)
hold off