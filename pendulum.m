% Cleanup step
clc;
clear;             
close all;

% Boundary values
y = [pi/4; 10];
t = 0;
u = 0;

% Sim parameters
dt = 0.05;
t_end = 3.0*pi;

% Physics here
g = -9.81;   % m/s^2
l = 1;      % m

% LQR  weigts
Q = [1 0; 0 1];
R = [1];

% Main loop
for i = 1 : t_end/dt
    [A B] = jacob(@RHS, y, t, u);
    
    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
    y=aa_rk45(@RHS, y, t, dt, u);
    
    % Plotting
    tp(i) = t;      % Time vector - for plotting
    up(i) = u;      % Control vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting
    refresh;
    % plot(tp, yp(1,:), 'r', tp, yp(2,:), 'b', tp, up, 'g', tp(i), up(i), 'xg', tp(i), yp(1,i), 'or', tp(i), yp(2,i), 'xb');
    % pause(0.001)
    % txt = sprintf('Symulacja dzialania wahadla matematycznego');
    % grid on;
    % title(txt);
  
    t=t+dt;
end

plot(tp, yp(1,:), '--r', tp, yp(2,:), '--b')
txt = sprintf('Symulacja dzialania wahadla matematycznego');
grid on;
title(txt);