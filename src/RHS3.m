% Right hand side function for point-in-space spacecraft
% y = [du dw dx dz]

function dy = RHS3(y, t, u)
R = 7080e3;  % m     - geostacionary orbit
g = 7.9536;  % m/s^2 - gravity at r
m = 2623;    % kg    - mass of the spacecraft
n = sqrt(g/R);

A = [0     0     n     -n^2  0     0     ;
     0     0     0     0     -n^2  0     ;
     -n    0     0     0     0     2*n^2 ;
     1     0     0     0     0     n     ; 
     0     1     0     0     0     0     ;
     0     0     1     -n    0     0    ];

B = [1     1     1     0     0     0    ]'/m;

dy = A*y + B*u;
end
