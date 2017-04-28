 function my_ode()
 options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
 [T,C]   = ode45(@test,[0 1],[0.021 0.0105 0.5],options);
 plot(T,C(:,1),'-',T,C(:,2),'-.',T,C(:,3),'.')
 pause(0.01)
 end

 function dc = test(t,c)
 k1    = 55.2;
 k2    = 30.2;
 dc    = zeros(3,1);
 dc(1) = -k1*c(1)^(1/2)*c(2) - k2*c(3)*c(1)^(1/2);
 dc(2) = -k1*c(2)*c(1)^(1/2);
 dc(3) =  k1*c(2)*c(1)^(1/2) - k2*c(3)*c(1)^(1/2);
 end