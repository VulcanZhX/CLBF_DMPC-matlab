f = @(x, u) -x + u;
[t, x_sol]=ode45(@(t, x)f(x, 1) , [0 1], 2);
x_sol_fixed = interp1(t, x_sol, (0:0.0001:1));
discreted_normed_integeral = 0.0001*sum(trace(x_sol_fixed'*x_sol_fixed));
disp(discreted_normed_integeral);