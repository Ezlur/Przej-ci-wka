% Cleanup step
clc;
clear;             
close all;

% Boundary values
#y = [1 -1 -6 50 -100 350]';
y = [0 0 0 0 0 0]';
t = 0;
u = [0 0 0]';

% Simulation parameters
dt = 1;
t_end = 600;
eps = 0.1;
deps = 0.01;
lol = 1;

% LQR weights
Q = 1*eye(6);
R = 1*eye(3);

% Main loop
for i = 1 : t_end/dt
    [A B] = jacob(@RHS2, y, t, u);
    [K P] = lqr_m(A, B, Q, R);
    u = -K*y;
%    u = [0 0 0]';
    y=aa_rk45(@RHS2, y, t, dt, u);
    
    if t==60
      y(4)=y(4)+0.015424; 
    end
    
    tp(i) = t;      % Time vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting
    up(:,i) = u;    % Control vector - for plotting

    if y(1)<deps && y(2)<deps && y(3)<deps && y(4)<eps && y(5)<eps && y(6)<eps && lol==1
      t
      lol=0;
    end

    t=t+dt;
end

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

print(h,'-dpng','-color','vib_plt.png')