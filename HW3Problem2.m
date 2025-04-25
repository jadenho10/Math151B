% adaptive_rk2_euler.m
clear; clc; close all;

%Define the ODE and its exact solution
% ODE: du/dx = -1000*(u - cos(x)) + sin(x)
f = @(x,u) -1000*(u - cos(x)) - sin(x);
% Exact solution: u(x) = cos(x)
u_exact = @(x) cos(x);

%Problem parameters
x0 = 0;       % initial x
xf = 10;      % final x
u0 = 1;       % initial condition: u(0) = 1
tol = 1e-6;   % tolerance for the adaptive method

% Initial step size
h = 0.002;

%Initialize storage arrays for adaptive solution and step sizes
% x_vals and u_vals will store the accepted (x, u) values
x_vals = x0;
u_vals = u0;
% h_vals will store the step sizes that were used (each corresponding to an accepted step)
h_vals = [];

% Initialization of integration variables
t = x0;
u_current = u0;

%Adaptive integration using embedded RK2 and Euler method
while t < xf
    % Adjust the step size if the next step would exceed xf
    if t + h > xf
        h = xf - t;
    end
    
    % Compute one Euler step (first-order method)
    u_euler = u_current + h * f(t, u_current);
    
    % Compute one RK2 (midpoint) step (second-order method)
    k1 = f(t, u_current);
    k2 = f(t + h/2, u_current + (h/2)*k1);
    u_rk2 = u_current + h * k2;
    
    % Estimate the local error as the absolute difference between the two methods
    err = abs(u_rk2 - u_euler);
    
    % Adaptive decision: if error is within tolerance, accept the step
    if err <= tol
        % Advance the solution
        t = t + h;
        u_current = u_rk2;
        % Store the accepted (t, u) and the step size used
        x_vals(end+1) = t; %#ok<*SAGROW>
        u_vals(end+1) = u_current;
        h_vals(end+1) = h; 
    end
    
    % Update step size:
    % When err == 0 we simply double the step, otherwise compute the scaling factor.
    if err == 0
        factor = 2;
    else
        factor = 0.9 * sqrt(tol/err);  % exponent 1/2 because the error ~ O(h^2)
    end
    % Optionally restrict the factor to avoid extreme changes
    factor = min(max(factor, 0.2), 5.0);
    
    % Update the step size for the next iteration
    h = h * factor;
end

%(Optional) Compute the exact solution at the adaptive grid points for comparison
u_exact_vals = u_exact(x_vals);

%Compute the relative L2 norm error:
% relative error = ||u_adaptive - u_exact||_2 / ||u_exact||_2
relL2 = norm(u_vals - u_exact_vals, 2) / norm(u_exact_vals, 2);
fprintf('Relative L2 norm error: %e\n', relL2);

fprintf('Number of steps taken: %d\n', length(x_vals) - 1);

%Plot the numerical solution versus the exact solution
figure;
plot(x_vals, u_exact_vals, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
hold on;
plot(x_vals, u_vals, 'r--', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Adaptive RK2/Euler');
xlabel('x');
ylabel('u(x)');
title('Numerical vs Exact Solution');
legend('Location', 'best');
yl = ylim;
ylim([-1.5, yl(2)]);
grid on;
hold off;

% Plot the adaptive step size versus x
figure;
plot(x_vals(1:end-1), h_vals, 'b-', 'LineWidth', 1.5);
xlabel('x');
ylabel('Step size, h');
title('Step Size (h) vs. x');
grid on;
