% Euler's Method for du/dx = (u^2 + u)/x, u(1) = -2, with exact solution
clc; clear; close all;

% Define function and exact solution
f = @(x, u) (u.^2 + u) ./ x;
exact = @(x) (2*x) ./ (1 - 2*x); 

x_start = 1;
x_end = 3;

% Step sizes to test
step_sizes = [0.5, 0.2, 0.1];
colors = {'r', 'g', 'b'};
legend_entries = {};

figure;
hold on;

% Loop over each step size
for k = 1:length(step_sizes)
    h = step_sizes(k);
    x = x_start:h:x_end;
    u = zeros(size(x));
    u(1) = -2;

    % Euler iteration
    for j = 1:(length(x) - 1)
        u(j+1) = h * f(x(j), u(j)) + u(j);
    end

    % Plot numerical solution
    plot(x, u, colors{k}, 'DisplayName', ['Euler, h = ', num2str(h)]);
    legend_entries{end+1} = ['Euler, h = ', num2str(h)];
end

% Plot exact solution
x_exact = linspace(x_start, x_end, 1000);
plot(x_exact, exact(x_exact), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exact');
legend_entries{end+1} = 'Exact';

legend('show');
xlabel('x');
ylabel('u(x)');
title('Euler''s Method vs Exact Solution');
grid on;
