% Cleanup step
clc;
clear;             
close all;

% Boundary values
y = [10 200 21 -48]';
t = 0;
u = [0 0 0 0]';

% Simulation parameters
dt = 0.1;
t_end = 10;

% LQR weights
Q = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
R = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

% Main loop
for i = 1 : t_end/dt
    
    [A B] = jacob(@RHS2, y, t, u);
    
    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
    y=aa_rk45(@RHS2, y, t, dt, u);
    
    % Plotting
    tp(i) = t;      % Time vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting
%    refresh;
%    plot(tp, yp(3,:), 'r', tp, yp(4,:), 'b', tp(i), yp(1,i), 'or', tp(i), yp(2,i), 'xb');
%    pause(0.001)
%    txt = sprintf('Symulacja dzialania punktowego satelity');
%    grid on;
%    title(txt);
  
    t=t+dt;
end

y

plot(tp, yp(4,:), '--r', tp, yp(3,:), '--b')
txt = sprintf('Symulacja dzialania punktowego satelity');
grid on;
title(txt);