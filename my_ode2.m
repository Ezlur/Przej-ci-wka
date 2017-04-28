function my_ode2()
clear  % Wyczyszczenie pamięci

t_end = 50*pi;
dt = 0.1;
t = 0;
%n = 2;
% Warunki poczatkowe - na razie niepotrzebne
y = [ 30; 1.0 ];  
%x=2;
u=0;
for i = 1 : t_end/dt
  y=aa_rk45(@fun, y, t, dt, u);
  tp(i) = t;      % Tablica czasow - do wykresu
  yp(:,i) = y;
  
  % Czesc merytoryczna - tu będzie wstawine całkowanie
  % x(1) = sin(t);
  % x(2) = cos(t);
  %x = aa_rk45(RHS,x,t,dt,u)
  % Tworzenie wykresow
  %tp(i) = t;      % Tablica czasow - do wykresu
  %yp(i) = x;  
  %yp(:,i) = x;    % Tablica wartosci funkcji - do wykresu
  refresh;
  plot( tp , yp(1,:) , 'r' , tp , yp(2,:) , 'b' , tp(i) , yp(1,i) , 'or' , tp(i) , yp(2,i) , 'xb' );
  pause(0.01)
  grid on
  txt = sprintf('Symulacja sinusa i cosinusa, t = %10.4f    y(1)=%10.5f  y(2)=%10.5f' , t , y(1) , y(2));
  title( txt );

  t = t + dt;
end

end

function dy = fun(y,t,u)
g=9.81; %m/s2
l=1; %m
dy=[y(2); -(g/l)*sin(pi*y(1)/180)];
end

function y = aa_rk45(RHS, x, t, dt, u)
  y = zeros(length(x), 1);

  y0 = feval( RHS , x , t, u );

  t1 = t + dt*0.25;
  vec = x + dt*(0.25)*y0;
  y1 = feval( RHS, vec , t1, u );

  t2 = t + dt*(3.0/8.0);
  vec = x + dt*( (3.0/32.0)*y0 + (9.0/32.0)*y1 );
  y2 = feval( RHS, vec , t2, u );

  t3 = t + dt*(12.0/13.0);
  vec = x + dt*( (1932.0/2197.0)*y0 + (-7200.0/2197.0)*y1 + (7296.0/2197.0)*y2 );
  y3 = feval( RHS, vec , t3, u );

  t4 = t + dt;
  vec = x + dt*( (439.0/216.0)*y0 + (-8.0)*y1 + (3680.0/513.0)*y2 + (-845.0/4104.0)*y3 );
  y4 = feval( RHS, vec , t4, u );

  t5 = t + dt*(1.0/2.0);
  vec = x + dt*( -(8.0/27.0)*y0 + (2.0)*y1 + (-3544.0/2565.0)*y2 + (1859.0/4104.0)*y3 + (-11.0/40.0)*y4 );
  y5 = feval( RHS, vec , t5, u );

  y = x + dt * ( (16.0/135.0)*y0 + (6656.0/12825.0)*y2 + (28561.0/56430.0)*y3 + (-9.0/50.0)*y4 + (2.0/55.0)*y5 );

end


