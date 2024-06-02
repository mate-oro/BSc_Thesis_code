% folder path for functions
addpath('methods/');

% rhs f
f = @(t,y) ([-0.04*y(1) + 10^4*y(2)*y(3);
           0.04*y(1) - 10^4*y(2)*y(3) - 3*10^7*y(2)^2;
                                        3*10^7*y(2)^2]);

% initial setup
a = 0;
b = 0.25;
y0 = [1; 0; 0];

% solving the IVP
options1 = odeset('RelTol', 1e-6,'AbsTol',1e-6,'Stats','on');
[T1, D1] = ode45(f,[0 b], y0, options1);
[T2, D2] = DOPRI5(f,0,b+0.02, y0, 1e-6);
[T3, I3] = AGauss2(f, 0, b, y0, 1e-6);
[T4, D4] = DOPRI5(f,0,b+0.02, y0, 2e-7);

% plotting the results, y(2)
figure();
subplot(3,1,[1 2]);
hold on;
plot(T1, D1(:,2),'-', 'Color','#00cd6c', 'MarkerFaceColor','#00cd6c','MarkerSize',3);
plot(T2, D2(2,:),'o-', 'Color','#009ade', 'MarkerFaceColor','#009ade','MarkerSize',5); 
plot(T4, D4(2,:),'o-', 'Color','#af58ba', 'MarkerFaceColor','#af58ba','MarkerSize',5);
plot(T3, I3(2,:),'o-', 'Color','#ffc61e', 'MarkerFaceColor','white','MarkerSize',7,'LineWidth', 1.5); 
ylabel('y(2)');
legend('MATLAB ode45', 'DOPRI','DOPRI pontos', 'Gauss2','northwest','FontSize',10);
ylim((1e-5)*[3.4 3.7]);
xlim([0 b]);

% number of steps for each method
disp("ode45 steps: " + max(size(find(T1 <= b))));
disp("DOPRI steps: " + max(size(find(T2 <= b))));
disp("DOPRI pontos steps: " + max(size(find(T4 <= b))));
disp("AGauss2 steps: " + max(size(find(T3 <= b))));

% calculating step sizes
H1 = T1(1:max(size(find(T1 <= b)))-1);
H2 = T2(1:max(size(find(T2 <= b)))-1);
H3 = T3(1:max(size(find(T3 <= b)))-1);
H4 = T4(1:max(size(find(T4 <= b)))-1);
for i=1:max(size(H1))
    H1(i) = T1(i+1) - T1(i);
end
for i=1:max(size(H2))
    H2(i) = T2(i+1) - T2(i);
end
for i=1:max(size(H3))
    H3(i) = T3(i+1) - T3(i);
end
for i=1:max(size(H4))
    H4(i) = T4(i+1) - T4(i);
end

% plotting the step size
subplot(3,1,3);
hold on;
plot(T1(1:max(size(find(T1 <= b)))-1),H1,'.-', 'Color','#00cd6c','LineWidth', 1.25);
plot(T4(1:max(size(find(T4 <= b)))-1),H4,'.-', 'Color','#af58ba','LineWidth', 1.25);
plot(T2(1:max(size(find(T2 <= b)))-1),H2,'.-', 'Color','#009ade','LineWidth', 1.25);
plot(T3(1:max(size(find(T3 <= b)))-1),H3,'.-', 'Color','#ffc61e','LineWidth', 1.25);
xlabel('t');
ylabel('h');
xlim([1e-8 b]);
set(gca,'Yscale', 'log');