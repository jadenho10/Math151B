% rk2_basic_solver_fixed.m
clear; clc; close all;

% Define the ODE as an anonymous function.
% Equation: du/dx = -1000*(u - cos(x)) + sin(x)
f = @(x, u) -1000*(u - cos(x)) - sin(x);

% Define the exact solution.
u_exact = @(x) cos(x);

% Problem parameters
x0 = 0;       % Starting x value
xf = 10;      % Ending x value
u0 = 1;       % Initial condition: u(0) = 1
h = 0.002;    % Fixed step size

% Determine the number of steps and create the x-grid.
N = floor((xf - x0) / h);
x = x0:h:xf;

% Preallocate the solution vector and initialize the initial condition.
u = zeros(size(x));
u(1) = u0;

% RK2 (Midpoint) method:
for j = 1:N
    k1 = f(x(j), u(j));                   % Slope at the current point.
    k2 = f(x(j) + h/2, u(j) + (h/2)*k1);    % Slope at the midpoint.
    u(j+1) = u(j) + h * k2;                % Update the solution.
end

% Compute the exact solution and calculate the relative L2 norm error.
u_ex = u_exact(x);
relL2 = norm(u - u_ex, 2) / norm(u_ex, 2);
fprintf('For h = %g, Relative L2 norm error = %e\n', h, relL2);

% Plot the exact solution and the RK2 numerical solution.
figure;
plot(x, u_ex, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
hold on;
plot(x, u, 'r--', 'LineWidth', 1.5, 'DisplayName', 'RK2 Numerical Solution');
xlabel('x');
ylabel('u(x)');
title(sprintf('RK2 (Midpoint) Numerical vs Exact Solution for h = %g', h));
legend('show', 'Location', 'best');
yl = ylim;
ylim([-1.5, yl(2)]);
grid on;
