% Numerical jacobian calculation

function [A B] = jacob(RHS, y, t, u)
delta = 1.0e-6;

for i = 1:size(y)(1)
  dy = zeros(size(y)(1), 1);
  dy(i) = delta;
  A(:,i) = (feval(RHS, y+dy, t, u) - feval(RHS, y, t, u))/delta;
endfor

for i = 1:size(u)(1)
  du = zeros(size(u)(1), 1);
  du(i) = delta;
  B(:,i) = (feval(RHS, y, t, u+du) - feval(RHS, y, t, u))/delta;
endfor

end