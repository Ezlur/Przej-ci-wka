% Numerical jacobian calculation

function M = jacob(RHS, y, t, u)

delta = 1.0e-6;

for i = 1:size(y)(1)
  dy = zeros(size(y)(1), 1);
  dy(i) = delta;
  M(:,i) = (feval(RHS, y+dy, t, u) - feval(RHS, y, t, u))/delta;
endfor

end