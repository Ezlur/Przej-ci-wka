clear  % Wyczyszczenie pamięci
close all

t_end = 4.0*pi;
dt = 0.1;
t = 0;

n = 2;
%x = [ 0.0 , 1.0 ];   % Warunki początkowe - na razie niepotrzebne

for i = 1 : t_end/dt

  % Część merytoryczna - tu będzie wstawine całkowanie
   x(1) = sin(t);
   x(2) = cos(t);
  %x=aa_rk45(RHS, x, t, dt, u)
  

  % Tworzenie wykresów
  tp(i) = t;      % Tablica czasów - do wykresu
  yp(:,i) = x;    % Tablica wartości funkcji - do wykresu
  % refresh;
  plot( tp , yp(1,:) , 'r' , tp , yp(2,:) , 'b' , tp(i) , yp(1,i) , 'or' , tp(i) , yp(2,i) , 'xb' );
  pause(0.01)
  txt = sprintf('Symulacja sinusa i cosinusa t =10.4f x(1)=10.5f x(2)=10.5f x(1) , x(2)');
  title(txt);

  t = t + dt;
end

