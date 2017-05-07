% Right hand side function for point-in-space spacecraft

function dy = RHS2(y, t, u)
g = 9.81;    % m/s^2 - gravity
R = 42160e3; % m     - geostacionary orbit
m = 1000;    % kg    - mass of the spacecraft
n = sqrt(g/R);

A = [0     n     -n^2  0    ;
     -n    0     0     2*n^2;
     1     0     0     n    ; 
     0     1     -n    0   ];

B = [1     1     1     1   ]'/m;

dy = A*y + u;
end
