% Cleanup step
clc;
clear;             
close all;

% Boundary values
y = [10 200 -100 21 -48 21]';
t = 0;
u = [0 0 0]';

% Simulation parameters
dt = 0.1;
t_end = 1000;

% LQR weights
Q = eye(6);
R = eye(3);

% Main loop
for i = 1 : t_end/dt

    [A B] = jacob(@RHS2, y, t, u);
    
    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
    y=aa_rk45(@RHS2, y, t, dt, u);
    
    tp(i) = t;      % Time vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting

    t=t+dt;
end

K
u

% Plotting
plot(tp, yp(4,:), '--g', tp, yp(5,:), '--r', tp, yp(6,:), '--b')
txt = sprintf('Symulacja dzialania punktowego satelity');
grid on;
legend('x', 'y', 'z');
xlabel('t [s]');
ylabel('odchylenie [m]');
title(txt);