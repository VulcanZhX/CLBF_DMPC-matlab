f = @(x, u) -eye(2)*x + eye(2)*u;
[t, x_sol]=ode45(@(t, x)f(x, [1 1]') , [0 1], [1 1]');
x_sol_fixed = interp1(t, x_sol, (0:1e-4:1-1e-4));
discreted_normed_integeral = 1e-4*sum(trace(x_sol_fixed'*x_sol_fixed));
disp(discreted_normed_integeral);