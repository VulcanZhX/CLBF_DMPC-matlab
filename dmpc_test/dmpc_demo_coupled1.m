%% 2-d cooperative overlapping decompositon iterative DMPC-demo without state constraints
clear
rng('default')
global sample_interval simulation_interval Nc Np Q R
% define LTI system:
A = [-1 0.3; -0.15 -1.5];
B = [0.6 0.1;0.2 0.8];
sys = @(x, u) A*x + B*u; % system function


% initial state and target state
x0_initial = [2.25 0.5]'; xs = [0 0]'; x0 = x0_initial;

% time parameters
sample_interval = 0.05; simulation_interval = 1e-4;
sim_steps = 1000; % number of simulation steps, overall time = sim_steps*sample_interval
max_iteration = 100; % number of iterations within one sample interval

% controller parameters
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = eye(2); % state cost
R = 1e-3*eye(2); % control cost

% initial input guess and bounds
U_initial = rand(Nc, 2);
U_lb = -2*ones(Nc, 2); % lower bound
U_ub = 2*ones(Nc, 2); % upper bound

% log the state and control input
x_log = []; u_log = [];

% barriers
x_barrier = [1 0]'; d_barrier = 0.5;

% rolling optimization
opt_options = optimoptions('fmincon', 'Display', 'off');
simulation_timer = 0;
for i_sim = 1:sim_steps
    tic; % begin current loop timer
    for i_p = 1:max_iteration
        % optimization problem for input 1
        U_host_init = U_initial(:, 1);
        U_adj_init = U_initial(:, 2);
        localcost = @(U)sub_cost_function(sys, U, U_adj_init, x0);
        nonlin_handle = @(U)nonlinear(sys, U, U_adj_init, x0, x_barrier, d_barrier);
        opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'nonlcon', nonlin_handle, 'options', opt_options);
        % opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
        %   'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'options', opt_options);
        [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input 2
        U_host_init = U_initial(:, 2);
        U_adj_init = U_initial(:, 1);
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 2), 'ub', U_ub(:, 2), 'nonlcon', nonlin_handle, 'options', opt_options);
        % opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
        %    'x0', U_host_init, 'lb', U_lb(:, 2), 'ub', U_ub(:, 2),  'options', opt_options);
        [U_opt2, J_opt2] = fmincon(opt_problem2);

        % convex combination of the two inputs
        weightcost_obj = @(gamma) sub_cost_function(sys, gamma(1)*U_opt1 + (1-gamma(1))*U_initial(:, 1),...
            gamma(2)*U_opt2 + (1-gamma(2))*U_initial(:, 2), x0); % gamma in [0, 1]

        % nonlinear constraint test
        nonlin_gamma_cons = @(gamma)nonlinear_gamma(sys, gamma, U_opt1, U_opt2, U_initial(:, 1), U_initial(:, 2), x0, x_barrier, d_barrier);
        opt_problem_convex_weight = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
            'x0', [0.5 0.5], 'lb', [0 0], 'ub', [1 1], 'nonlcon', nonlin_gamma_cons, 'options', opt_options);
        % opt_problem_convex_weight1 = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
        %    'x0', [0.5 0.5], 'lb', [0 0], 'ub', [1 1], 'options', opt_options);
        [gamma_opt, J_opt] = fmincon(opt_problem_convex_weight);
        U_opt1_new = gamma_opt(1)*U_opt1 + (1-gamma_opt(1))*U_initial(:, 1);
        U_opt2_new = gamma_opt(2)*U_opt2 + (1-gamma_opt(2))*U_initial(:, 2);

        % iteration terminal condition
        if norm(U_opt1_new - U_initial(:, 1)) < 1e-2 && norm(U_opt2_new - U_initial(:, 2)) < 1e-2
            break;
        end
        % update control input
        U_initial = [U_opt1_new U_opt2_new];
    end
    % update state
    [t, x_full_sol] = ode45(@(t, x) sys(x, U_initial(1, :)'), [0 sample_interval], x0);
    x0_new = x_full_sol(end, :)'; % new state
    x_log = [x_log x0_new]; % log the state
    u_log = [u_log U_initial(1, :)']; % log the control input
    x0 = x0_new; % update state
    % display elasped time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f\n', i_sim, toc, norm(x0 - xs));
    % time forwarding terminal condition
    if norm(x0 - xs) < 0.5
        sample_interval = 0.005;
        R = 1e-2*eye(2);
    end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0 - xs) < 0.02
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end

% simulation time
fprintf('Total simulation time: %.4f seconds\n', simulation_timer);
% plot state
figure
plot(x_log(1,:), x_log(2,:), 'LineWidth', 1, 'Color', 'b')
hold on
fimplicit(@(x1,x2) (x1-x_barrier(1))^2+(x2-x_barrier(2))^2-d_barrier^2, "--")
axis([-1 4 -1 4])
grid on
xlabel('x1')
ylabel('x2')
title('State evolution')


function J = sub_cost_function(subsys, u_host, u_adj, x0)
% sub_cost_function: cost function for the subsystem
global sample_interval simulation_interval Nc Np Q R

% remember to delete this
x0_host_s = [0 0]'; % target state

J = 0; % initialize cost

u_overall = [u_host u_adj];

for i = 1:Np
    if i <= Nc
        u_i = u_overall(i, :)';
    else
        u_i = u_overall(Nc, :)';  % Nc <= Np
    end
    [t, sub_xhost_sol] = ode45(@(t, x) subsys(x, u_i), [0 sample_interval], x0);
    sub_xhost_fixed_interval = interp1(t, sub_xhost_sol, (0:simulation_interval:sample_interval-simulation_interval));

    % calculate cost in one sample interval
    delta_subx_array = sub_xhost_fixed_interval - repmat(x0_host_s', size(sub_xhost_fixed_interval, 1), 1); % X - Xs
    sub_weighted_norm = trace(delta_subx_array * Q * delta_subx_array');
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i' * R * u_i * simulation_interval;
    x0 = sub_xhost_fixed_interval(end, :)';
end
end


function [cineq, ceq] = nonlinear(subsys, u_host, u_adj, x0, x_barrier, d_barrier)
% suppose ||x0-x_barrier||>=d_barrier
global sample_interval

Np_nonlin = 1;
cineq = zeros(Np_nonlin, 1); ceq = [];
u_overall = [u_host u_adj];
u_in = u_overall(1, :)';

[t, x_sol] = ode45(@(t,x) subsys(x, u_in), [0 Np_nonlin*sample_interval], x0);
t_sample = sample_interval:sample_interval:Np_nonlin*sample_interval;
x_sol_sample = interp1(t, x_sol, t_sample);
cineq(1, :) = d_barrier - norm(x_sol_sample(1,:)'-x_barrier);
end

function [cineq, ceq] = nonlinear_gamma(subsys, gamma, u_host, u_adj, u_host_pre, u_adj_pre, x0, x_barrier, d_barrier)
% suppose ||x0-x_barrier||>=d_barrier
global sample_interval
% gamma = [gamma(1) gamma(2)]; gamma(1) is the weight of u_host, gamma(2) is the weight of u_adj
Np_nonlin = 1;
cineq = zeros(Np_nonlin, 1); ceq = [];
u_host = gamma(1)*u_host + (1-gamma(1))*u_host_pre;
u_adj = gamma(2)*u_adj + (1-gamma(2))*u_adj_pre;
u_overall = [u_host u_adj];

u_in = u_overall(1, :)';

[t, x_sol] = ode45(@(t,x) subsys(x, u_in), [0 Np_nonlin*sample_interval], x0);
t_sample = sample_interval:sample_interval:Np_nonlin*sample_interval;
x_sol_sample = interp1(t, x_sol, t_sample);
cineq(1, :) = d_barrier - norm(x_sol_sample(1,:)'-x_barrier);
end