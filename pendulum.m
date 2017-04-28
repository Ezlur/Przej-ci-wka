% Cleanup step
clc;
clear;             
close all;

% Boundary values
y = [pi/4; 0];
t = 0;
u = 0;

% Sim parameters
dt = 0.1;
t_end = 10.0*pi;

% Physics here
g = 9.81;   % m/s^2
l = 1;      % m

% LQR  weigts
Q = [1 0; 0 1];
R = [1];

% Main loop
for i = 1 : t_end/dt
    A = [0 1; -g*(cos(y(1)))/l 0];
    B = [0; 1];
    
    [K P] = lqr_m(A, B, Q, R)
    u = -K*y;
    
    y=aa_rk45(@RHS, y, t, dt, u);
    
    % Plotting
    tp(i) = t;      % Time vector - for plotting
    up(i) = u;      % Control vector - for plotting
    yp(:,i) = y;    % Function matrix - for plotting
    refresh;
    plot(tp, yp(1,:), 'r', tp, yp(2,:), 'b', tp, up, 'g', tp(i), up(i), 'xg', tp(i), yp(1,i), 'or', tp(i), yp(2,i), 'xb');
    pause(0.001)
    txt = sprintf('Symulacja dzia³ania wachadla matematycznego');
    grid on;
    title(txt);
  
    t=t+dt;
end

% Right hand side function
function dy = RHS(y,t,u)

g = 9.81;   % m/s^2
l = 1;      % m
dy = [y(2); -(g/l)*sin(pi*y(1)/180) + u];
end

% Numerical integration
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

% Linear Quadratic Regulator
function [k,s]=lqr_m(a,b,q,r,nn)

%LQR	Linear quadratic regulator design for continuous-time systems.
%	[K,S] = LQR(A,B,Q,R)  calculates the optimal feedback gain matrix K
%	such that the feedback law  u = -Kx  minimizes the cost function:
%
%		J = Integral {x'Qx + u'Ru} dt
%
%	subject to the constraint equation: 
%		.
%		x = Ax + Bu 
%                
%	Also returned is S, the steady-state solution to the associated 
%	algebraic Riccati equation:
%				  -1
%		0 = SA + A'S - SBR  B'S + Q
%
%	[K,S] = LQR(A,B,Q,R,N) includes the cross-term N that relates
%	u to x in the cost function.

%	J.N. Little 4-21-85
%	Revised 8-27-86 JNL
%	Copyright (c) 1985, 1986 by the MathWorks, Inc.

%error(nargchk(4,5,nargin));
%error(abcdchk(a,b));

[m,n] = size(a);
[mb,nb] = size(b);
[mq,nq] = size(q);
if (m ~= mq) || (n ~= nq) 
	error('A and Q must be the same size')
end
[mr,nr] = size(r);
if (mr ~= nr) || (nb ~= mr)
	error('B and R must be consistent')
end

if nargin == 5
	[mn,nnn] = size(nn);
	if (mn ~= m) || (nnn ~= nr)
		error('N must be consistent with Q and R')
	end
	% Add cross term
	q = q - nn/r*nn';
	a = a - b/r*nn';
else
	nn = zeros(m,nb);
end

% Check if q is positive semi-definite and symmetric
if any(eig(q) < 0) || (norm(q'-q,1)/norm(q,1) > eps)
	error('Q must be symmetric and positive semi-definite')
end
% Check if r is positive definite and symmetric
if any(eig(r) <= 0) || (norm(r'-r,1)/norm(r,1) > eps)
	error('R must be symmetric and positive definite')
end

% Start eigenvector decomposition by finding eigenvectors of Hamiltonian:
[v,d] = eig([a b/r*b';q, -a']);
d = diag(d);
[d,index] = sort(real(d));	 % sort on real part of eigenvalues

%if (~( (d(n)<0) && (d(n+1)>0) ))
if (~( (d(n)<1.e-15) && (d(n+1)>-1.e-15) ))
   printf('Can''t order eigenvalues, (A,B) may be uncontrollable. Checking rank C = [B A*B ... A^(n-1)*B ]\n')
   C = zeros(m,n*nb);
   c = b;
   C(:,1:nb) = c;
   for i = 1 : n-1 ,
      c = a*c;
      C(:,i*nb+1:i*nb+nb) = c;
   end
   %C
   rank_C = rank(C);
   if( rank_C < n )
     rank_C
     error('rank_C < n , (A,B) are uncontrollable.')
   else
     error('rank_C = n (OK), but there is something wrong with ordering - check it out!')
   end
end

chi = v(1:n,index(1:n));	 % select vectors with negative eigenvalues
lambda = v((n+1):(2*n),index(1:n));
s = -real(lambda/chi);
k = r\(nn'+b'*s);

end