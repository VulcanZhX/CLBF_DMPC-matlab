clear;
rng("default"); % For reproducibility

% system definition: x' = Ax + Bu (4D system with 4D input)
A = [-1 0.3 0 0;
    0.15 -1.5 0.1 1;
    0 0 -0.75 0.01;
    0 0 0 -1.5];

B = blkdiag([0.6 0.1; 0.2 0.8], [0.9 0.1; 0.2 0.8]);

% global func1 func2
global sample_interval simulation_interval Nc Np

func = @(x, u) A*x + B*u; % system function
func1 = @(x1, x2, u1, u2) A(1:2, :)*[x1; x2]  + B(1:2,1:2)*u1;
func2 = @(x1, x2, u1, u2) A(3:4, :)*[x1; x2]  + B(3:4,3:4)*u2;

% initial state and target state
x0_initial = [1.5 1 -2 -2]'; xs = [0 0 0 0]'; x0 = x0_initial;

% time parameters
sample_interval = 0.02; simulation_interval = 1e-4;
sim_steps = 300; % number of simulation steps, overall time = sim_steps*sample_interval

% controller parameters
Nc = 2; % controller horizon
Np = 6; % prediction horizon
Q = eye(4); % state cost
R = 1e-3*eye(4); % control cost

% initial input guess and bounds
U_initial = rand(Nc, 4);
U_lb = -2*ones(Nc, 4); % lower bound
U_ub = 2*ones(Nc, 4); % upper bound

% log the state and control input
x_log = []; u_log = [];

% rolling optimization
opt_options = optimoptions('fmincon', 'Display', 'off');
simulation_timer = 0;
for i_sim = 1:sim_steps
    tic; % begin current loop timer
    % optimization problem for system 1
    sub_x0_struct = struct('xhost', x0(1:2), 'xhost_s', xs(1:2), 'xadj', x0(3:4), 'xadj_s', xs(3:4));
    U_host_init = U_initial(:, 1:2);
    U_adj_init = U_initial(:, 3:4);
    % function test
    localcost1 = @(U)sub_cost_function(func1, U, U_adj_init, sub_x0_struct, Q(1:2, 1:2), R(1:2, 1:2));
    opt_problem1 = createOptimProblem('fmincon', 'objective', localcost1, ...
        'x0', U_host_init, 'lb', U_lb(:, 1:2), 'ub', U_ub(:, 1:2), 'options', opt_options);
    [U_opt1, J_opt1] = fmincon(opt_problem1);

    % optimization problem for system 2
    sub_x0_struct = struct('xhost', x0(3:4), 'xhost_s', xs(3:4), 'xadj', x0(1:2), 'xadj_s', xs(1:2));
    U_host_init = U_initial(:, 3:4);
    U_adj_init = U_initial(:, 1:2);
    localcost2 = @(U)sub_cost_function(func2, U, U_adj_init, sub_x0_struct, Q(3:4, 3:4), R(3:4, 3:4));
    opt_problem2 = createOptimProblem('fmincon', 'objective', localcost2, ...
        'x0', U_host_init, 'lb', U_lb(:, 3:4), 'ub', U_ub(:, 3:4), 'options', opt_options);
    [U_opt2, J_opt2] = fmincon(opt_problem2);

    % update control input and state
    U_opt_now = [U_opt1(1, :) U_opt2(1, :)]';
    [t, x_global_sol] = ode45(@(t, x) func(x, U_opt_now), [0 sample_interval], x0);
    x0_new = x_global_sol(end, :)'; % new state
    x_log = [x_log x0_new]; % log the state
    u_log = [u_log U_opt_now]; % log the control input
    x0 = x0_new; % update state
    U_initial = [U_opt1 U_opt2];% update initial guess
    % display elasped time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f\n', i_sim, toc, norm(x0 - xs));
    % terminal condition
    if norm(x0 - xs) < 1
        sample_interval = 0.005;
    end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0 - xs) < 0.075
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end

% simulation time
fprintf('Total simulation time: %.4f seconds\n', simulation_timer);

% plot state
figure
subplot(2, 1, 1);
plot(x_log(1, :), x_log(2, :), '-*')
axis equal; grid on; axis([-3 3 -3 3])
title('State evolution');
xlabel('State x1');
ylabel('State x2');

subplot(2, 1, 2);
plot(x_log(3, :), x_log(4, :), '-*')
axis equal; grid on; axis([-3 3 -3 3])
title('State evolution');
xlabel('State x3');
ylabel('State x4');



function J = sub_cost_function(subsys, u_host, u_adj, sub_x0_struct, sub_Q, sub_R)
% sub_cost_function: cost function for the subsystem
global sample_interval simulation_interval Nc Np

J = 0; % initialize cost

% unpack the state and input variables
x0_host = sub_x0_struct.xhost;
x0_host_s = sub_x0_struct.xhost_s;
x0_adj = sub_x0_struct.xadj;

for i = 1:Np
    if i <= Nc
        u_i = u_host(i, :)';
        u_j = u_adj(i, :)';
    else
        u_i = u_host(Nc, :)';  % Nc <= Np
        u_j = u_adj(Nc, :)';
    end
    [t, sub_xhost_sol] = ode45(@(t, x) subsys(x, x0_adj, u_i, u_j), [0 sample_interval], x0_host);
    sub_xhost_fixed_interval = interp1(t, sub_xhost_sol, (0:simulation_interval:sample_interval-simulation_interval));

    % calculate cost in one sample interval
    delta_subx_array = sub_xhost_fixed_interval - repmat(x0_host_s', size(sub_xhost_fixed_interval, 1), 1); % X - Xs
    sub_weighted_norm = trace(delta_subx_array * sub_Q * delta_subx_array');
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i' * sub_R * u_i * simulation_interval;
    x0_host = sub_xhost_fixed_interval(end, :)';
end
end


function [cineq, ceq] = nonlinear(U, x0, Np, x_barrier, d_barrier)
% suppose ||x0-x_barrier||>=d_barrier
global func
global sample_interval

cineq = zeros(Np, 1); ceq = [];
U = U';
syms t;
uf(t) = [0; 0]*t;
for i = 1:Np
    if i < size(U, 2)
        uf(t) = uf(t) + U(:, i)*(heaviside(t-(i-1)*sample_interval) - heaviside(t-i*sample_interval));
    else
        uf(t) = uf(t) + U(:, size(U, 2))*(heaviside(t-(size(U, 2)-1)*sample_interval));
        break;
    end
end
uf_handle = matlabFunction(uf);

[t, x_sol] = ode45(@(t,x) func(x, uf_handle(t)), [0 Np*sample_interval], x0);
t_sample = sample_interval:sample_interval:Np*sample_interval;
x_sol_sample = interp1(t, x_sol, t_sample);
for i_sample = 1:Np
    cineq(i_sample, :) = d_barrier - norm(x_sol_sample(i_sample,:)'-x_barrier);
end
end