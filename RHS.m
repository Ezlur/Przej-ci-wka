% Right hand side function
function dy = RHS(y, t, u)
g = 9.81;   % m/s^2
l = 1;      % m
dy = [y(2); -(g/l)*sin(y(1)) + u];
end
