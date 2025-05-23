clear;
u_discrete = [2 0.9 2; 1.1 2.8 1];
x0 = [3 3]';
x_barrier = [1.5 1.5]';
d_barrier = 0.5;
syms t;
uf(t) = [0; 0]*t;
for i = 1:length(u_discrete)
    uf(t) = uf(t) + u_discrete(:, i)*(heaviside(t-(i-1)) - heaviside(t-i));
end

% uf2 = matlabFunction(uf);
% uf2(2.5)
% f = @(x, u) -x + u;
% [t, x_sol]=ode45(@(t, x)f(x, uf2(t)) , [0 length(u_discrete)], [3 3]');
% 
% figure;
% subplot(2,1,1)
% plot(t, x_sol);
% subplot(2,1,2)
% t = 0:0.001:length(u_discrete);
% plot(t, uf2(t));

global func
func = @(x, u) -x + u;
U = [1 0;0 1;1 0];
[cineq, ceq] = nonlinear(U, x0, 5, x_barrier, d_barrier);

function [cineq, ceq] = nonlinear(U, x0, Np, x_barrier, d_barrier)
    % suppose ||x0-x_barrier||>=d_barrier
    global func

    cineq = zeros(Np, 1); ceq = [];
    U = U';
    syms t;
    uf(t) = [0; 0]*t;
    for i = 1:Np
        if i < size(U, 2)
            uf(t) = uf(t) + U(:, i)*(heaviside(t-(i-1)) - heaviside(t-i));
        else
            uf(t) = uf(t) + U(:, size(U, 2))*(heaviside(t-(size(U, 2)-1)));
            break
        end
    end
    uf_handle = matlabFunction(uf);
    [t, x_sol] = ode45(@(t,x) func(x, uf_handle(t)), [0 Np], x0);
    t_sample = 1:Np;
    x_sol_sample = interp1(t, x_sol, t_sample);
    for i_sample = 1:Np
        cineq(i_sample, :) = d_barrier - norm(x_sol_sample(i_sample,:)'-x_barrier);
    end
end