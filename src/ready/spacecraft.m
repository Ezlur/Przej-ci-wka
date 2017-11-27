% Cleanup step
clc;
clear;             
close all;

% Boundary values
% y = [du dv dw dx dy dz]
y = [1 -1 -6 50 -100 350]';
t = 0;
u = [0 0 0]';

% Simulation parameters
dt = 1;       % simulation time step
t_end = 600;  % simulation time

% LQR weights
Q = 1*eye(6);
R = 1*eye(3);

% Main loop
for i = 1 : t_end/dt
    [A B] = jacob(@RHS, y, t, u);
    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
    y=aa_rk45(@RHS, y, t, dt, u);
    
    tp(i) = t;      % Time vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting
    up(:,i) = u;    % Control vector - for plotting

    t=t+dt;
    
    
% Plotting
    h=figure(1);

    subplot(2,1,1)
    plot(tp, yp(4,:), 'g', tp, yp(5,:), 'r', tp, yp(6,:), 'b')
    grid on;
    legend('x', 'y', 'z');
    xlabel('t [s]');
    ylabel('odchylenie [m]');
    
    subplot(2,1,2)
    plot(tp, up(1,:), 'g', tp, up(2,:), 'r', tp, up(3,:), 'b')
    grid on;
    legend('du', 'dv', 'dw');
    xlabel('t [s]');
    ylabel('przyspieszenie [m/s2]');
    refresh;
end

% Final plot
h=figure(1);

subplot(2,1,1)
plot(tp, yp(4,:), 'g', tp, yp(5,:), 'r', tp, yp(6,:), 'b')
grid on;
legend('x', 'y', 'z');
xlabel('t [s]');
ylabel('odchylenie [m]');

subplot(2,1,2)
plot(tp, up(1,:), 'g', tp, up(2,:), 'r', tp, up(3,:), 'b')
grid on;
legend('du', 'dv', 'dw');
xlabel('t [s]');
ylabel('przyspieszenie [m/s2]');