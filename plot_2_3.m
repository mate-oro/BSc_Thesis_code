% folder path for functions
addpath('methods/');

% rhs f
f = @(t,y) ([1 + y(1)^2*y(2) - 4*y(1);
           3*y(1) - y(1)^2*y(2)]);

% initial setup
a = 0;
b = 27;
y0 = [0; 0];

% solving the IVP
[T2, D2] = DOPRI5(f,0,b, y0, 1e-5);

% plotting the results
figure();
hold on;
plot(T2, 0*T2,'|', 'Color','#ffc61e','MarkerSize', 12, 'LineWidth', 1.2);
plot(T2, D2(1,:),'o-', 'Color','#00cd6c', 'MarkerFaceColor','white','MarkerSize',4,'LineWidth', 1.2);
plot(T2, D2(2,:),'o-', 'Color','#af58ba', 'MarkerFaceColor','white','MarkerSize',4,'LineWidth', 1.2);
xlim([0 b]);
pbaspect([3 2 1]);
legend('lépések','x komponens', 'y komponens','Location','northwest','FontSize',10);