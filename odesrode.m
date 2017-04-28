clear  % Wyczyszczenie pamięci

t_end = 4*pi;
dt = 0.1;
t = 0;
n = 2;
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

function dy = fun(y,t,u)
g=9.81; %m/s2
l=1; %m
dy=[y(2); -(g/l)*sin(pi*y(1)/180)];
end