% Parameters
sigma = 10;
rho   = 28;
beta  = 8 / 3;

% Right-hand side functions
g1 = @(t, x, y, z) sigma * (y - x);
g2 = @(t, x, y, z, h) x * (rho - z) - y;
g3 = @(t, x, y, z, h) x * y - beta * z;

% Time settings
t0 = 0;
t_final = 100;
h = 0.01;
t = t0:h:t_final;
N = length(t);

% Preallocate arrays for x, y, z
x = zeros(1, N);
y = zeros(1, N);
z = zeros(1, N);

% Initial conditions
x(1) = 1;
y(1) = 1;
z(1) = 1;

% RK2 integration loop
for j = 1:N-1

    % Stage 1
    k1 = g1(t(j), x(j), y(j), z(j));
    l1 = g2(t(j), x(j), y(j), z(j));
    m1 = g3(t(j), x(j), y(j), z(j));

    % Stage 2
    k2 = g1(t(j)+ h/2, x(j) + h/2 * k1, y(j) + h/2 * l1, z(j) + h/2 * m1);
    l2 = g2(t(j) + h/2, x(j) + h/2 * k1, y(j) + h/2 * l1, z(j) + h/2 * m1);
    m2 = g3(t(j) + h/2, x(j) + h/2 * k1, y(j) + h/2 * l1, z(j) + h/2 * m1);

    % Update
    x(j+1) = x(j) + h * k2;
    y(j+1) = y(j) + h * l2;
    z(j+1) = z(j) + h * m2;
end

% Plot 3D trajectory
figure;
plot3(x, y, z, 'b');
xlabel('x'); ylabel('y'); zlabel('z');
title('Lorenz System');
grid on;

% Plot time series
figure;
subplot(3,1,1);
plot(t, x, 'r');
xlabel('Time t'); ylabel('x(t)');
title('x(t) vs. t');

subplot(3,1,2);
plot(t, y, 'g');
xlabel('Time t'); ylabel('y(t)');
title('y(t) vs. t');

subplot(3,1,3);
plot(t, z, 'b');
xlabel('Time t'); ylabel('z(t)');
title('z(t) vs. t');

% Plot z vs x
figure;
plot(x, z, 'k');
xlabel('x'); ylabel('z');
title('z vs. x');
grid on;
