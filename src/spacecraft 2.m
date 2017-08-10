% Cleanup step
clc;
clear;             
close all;

% Boundary values
y = [10 200 -100 21 -48 21]';
t = 0;
u = [0 0 0 0 0 0]';

% Simulation parameters
dt = 0.1;
t_end = 10;

% LQR weights
Q = eye(6);
R = eye(1);

r = 7080e3;  % m     - geostacionary orbit
g = 7.9536;  % m/s^2 - gravity at r
m = 2623;    % kg    - mass of the spacecraft
n = sqrt(g/r);

A = [0     0     n     -n^2  0     0     ;
     0     0     0     0     -n^2  0     ;
     -n    0     0     0     0     2*n^2 ;
     1     0     0     0     0     n     ; 
     0     1     0     0     0     0     ;
     0     0     1     -n    0     0    ];
     
B = [1     1     1     0     0     0    ]'/m;

% Main loop
for i = 1 : t_end/dt

    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
    y=aa_rk45(@RHS3, y, t, dt, u);
    
    tp(i) = t;      % Time vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting

    t=t+dt;
end

% Plotting
plot(tp, yp(4,:), '--g', tp, yp(5,:), '--r', tp, yp(6,:), '--b')
txt = sprintf('Symulacja dzialania punktowego satelity');
grid on;
legend('x', 'y', 'z');
xlabel('t [s]');
ylabel('odchylenie [m]');
title(txt);