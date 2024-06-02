% folder path for functions
addpath('methods/');

% initial setup
rows = 3;
res = 0.005;
[real_part, imag_part] = meshgrid(-rows:res:rows, -rows:res:rows);
complex_points = complex(real_part, -imag_part);

% parameter for theta-method ar LENM2
a = 1/2;

% stability functions
Rt = @(z) abs((1 + (1-a)*z)/(1 - a*z));
R1 = @(z) abs((1 + z));
R2 = @(z) abs((1 + z + z^2/2));
R3 = @(z) abs((1 + z + z^2/2 + z^3/6));
R4 = @(z) abs((1 + z + z^2/2 + z^3/6 + z^4/24));
R_LENM2 = @(z) abs((2 + (2 - 2*a)*z)/(2 - 2*a*z + (2*a - 1)*z^2));

% calculating the stability function
Z = arrayfun(R4, complex_points);

%plotting the region
figure();
imagesc([-rows,rows],[rows,-rows], min(1,Z), [0 1])
hold on
contourf(real_part,imag_part, Z, [1 1],'b')
colorbar
colormap("hot")
grid
ax = gca;
ax.XTick = -rows:1:rows;
ax.YTick = -rows:1:rows;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off
