clear
rng(12345)

global m1 m2 m3 k0 hd kc
k0 = 1.1; kc = 0.25; hd = 0.3;
m1 = 1.5; m2 = 2; m3 = 1;

% system definition
func = @(x, u)[
    m1*x(2);
    -k0*x(1)*exp(-x(1))-kc*(x(1)-x(3))-hd*x(2)+u(1);
    m2*x(4);
    -k0*x(3)*exp(-x(3))-kc*(x(3)-x(1))-kc*(x(3)-x(5))-hd*x(4)+u(2);
    m3*x(6);
    -k0*x(5)*exp(-x(5))-kc*(x(5)-x(3))-hd*x(6)+u(3);
    ];

% initial state and target state
x0_initial = [0.4 0 0.4 0 0.5 0];
xs = [0 0 0 0 0 0]'; x0 = x0_initial;
% time parameters
global sample_interval simulation_interval
sample_interval = 0.01; simulation_interval = 1e-4;
sim_steps = 10000; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 10; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = diag([1 0.05 1 0.05 1 0.05]); % state cost
Ru = 1e-1* diag([0.1 1 0.1]); % control cost

% initial input guess and bounds
U_initial = 2 * rand(Nc, 3) - 1; % initial guess
U_lb = -5*ones(Nc, 3); % lower bound
U_ub = 5*ones(Nc, 3); % upper bound
% log the state and control input
x_log = []; u_log = [];

% barriers
x_barrier = [0.1 0.25]'; % barrier point
d_barrier = 0.1; % distance to barrier

% lyap_formulation


% terminal constriant cost function B(x) = d_barrier - ||x-x_barrier||
global B_terminal_cost
B_terminal_cost = @(x) max(5*(d_barrier - norm(x([3, 5])-x_barrier)), -0.2);

% rolling optimization
opt_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','interior-point','EnableFeasibilityMode',true);
simulation_timer = 0;
for i_sim = 1:sim_steps
    tic;
    for i_p = 1:max_p_iteration
        % optimization problem for input 1
        U_host_init = U_initial(:, 1);
        U_adj_init = U_initial(:, 2:3);
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 1);
        nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, x_barrier, d_barrier, Np, 1);
        opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'nonlcon', nonlin_handle, 'options', opt_options);
        [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input 2
        U_host_init = U_initial(:, 2);
        U_adj_init = [U_initial(:, 1) U_initial(:, 3)];
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 2);
        nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, x_barrier, d_barrier, Np, 2);
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 2), 'ub', U_ub(:, 2), 'nonlcon', nonlin_handle,  'options', opt_options);
        [U_opt2, J_opt2] = fmincon(opt_problem2);

        % optimization problem for input 3
        U_host_init = U_initial(:, 3);
        U_adj_init = [U_initial(:, 1) U_initial(:, 2)];
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 3);
        nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, x_barrier, d_barrier, Np, 3);
        opt_problem3 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 3), 'ub', U_ub(:, 3), 'nonlcon', nonlin_handle,  'options', opt_options);
        [U_opt3, J_opt3] = fmincon(opt_problem3);

        % convex combination of the three optimization problems
        weightcost_obj = @(gamma) sub_cost_function(func, gamma(1)*U_opt1 + (1-gamma(1))*U_initial(:, 1),...
            [gamma(2)*U_opt2 + (1-gamma(2))*U_initial(:, 2), gamma(3)*U_opt3 + (1-gamma(3))*U_initial(:, 3)], x0, xs, 1); % gamma in [0, 1]
        nonlin_handle_comb = @(gamma) nonlinear_gamma_horizon(func, gamma, U_opt1, [U_opt2 U_opt3], ...
        U_initial(:, 1), [U_initial(:, 2) U_initial(:, 3)], x0, x_barrier, d_barrier, Np, 1);

        % optimization problem for the convex combination
        opt_problem_comb = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
            'x0', [0.5; 0.5; 0.5], 'lb', [0; 0; 0], 'ub', [1; 1; 1], 'nonlcon', nonlin_handle_comb,  'options', opt_options);
        [gamma_opt, J_opt_comb] = fmincon(opt_problem_comb);
        % update control input
        U_opt1_new = gamma_opt(1)*U_opt1 + (1-gamma_opt(1))*U_initial(:, 1);
        U_opt2_new = gamma_opt(2)*U_opt2 + (1-gamma_opt(2))*U_initial(:, 2);
        U_opt3_new = gamma_opt(3)*U_opt3 + (1-gamma_opt(3))*U_initial(:, 3);

        % iteration terminal condition
        if norm([U_opt1_new U_opt2_new U_opt3_new] - U_initial) < 0.1
            break;
        end
        U_initial = [U_opt1_new U_opt2_new U_opt3_new]; % update initial guess
    end
    %update state
    [~, x_global_sol] = ode45(@(t, x) func(x, U_initial(1, :)'), [0 sample_interval], x0);
    x0_new = x_global_sol(end, :)'; % new state
    x_log = [x_log x0_new]; % log the state
    u_log = [u_log U_initial(1, :)']; % log the control input'
    x0 = x0_new; % update state
    % display elasped time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f, J_now: %.4f\n', i_sim, toc, norm(x0 - xs), J_opt_comb);
    % terminal condition
    % if norm(x0 - xs) < 0.25
    %     sample_interval = 0.005;
    %     Ru = 2e-3*eye(3);
    % else
    %     sample_interval = 0.01;
    %     Ru = 1e-5*eye(3);
    % end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0 - xs) < 0.02
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end
% simulation time
fprintf('Total simulation time: %.4f seconds\n', simulation_timer);
% plot state
figure;
grid on;
xlabel('Time (s)')
ylabel('State x')
title('State evolution')
subplot(2, 1, 1);
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(1, :), 'b-')
hold on
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(3, :), 'r-')
hold on
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(5, :), 'g-')
legend('r1', 'r2', 'r3')

subplot(2, 1, 2);
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(2, :), 'b-')
hold on
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(4, :), 'r-')
hold on
plot(0:sample_interval:sample_interval*(size(x_log, 2)-1), x_log(6, :), 'g-')
legend('v1', 'v2', 'v3')

figure;
grid on;
plot(0:sample_interval:sample_interval*(size(u_log, 2)-1), u_log(1, :), 'b-')
hold on
plot(0:sample_interval:sample_interval*(size(u_log, 2)-1), u_log(2, :), 'r-')
hold on
plot(0:sample_interval:sample_interval*(size(u_log, 2)-1), u_log(3, :), 'g-')
xlabel('Time (s)')
ylabel('Control input u')
title('Control input evolution')
legend('u1', 'u2', 'u3')

figure;
plot(x_log(3, :), x_log(5, :), 'b-')
hold on
fimplicit(@(x1,x2) (x1-x_barrier(1))^2+(x2-x_barrier(2))^2-0.95*d_barrier^2, "--")
axis equal

%% SUBSYSTEM COST FUNCTION
function J = sub_cost_function(subsys, u_host, u_adj, x0, xs, ORDER)
% sub_cost_function: cost function for the subsystem
global sample_interval simulation_interval Nc Np Q Ru
global B_terminal_cost
% remember to delete this
x0_host_s = xs; % target state

J = 0; % initialize cost
P_terminal = 2*Q;
if ORDER == 1
    u_overall = [u_host u_adj];
elseif ORDER == 2
    u_overall = [u_adj(:,1) u_host u_adj(:,2)];
else
    u_overall = [u_adj u_host];
end

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
    J = J + sub_discreted_normed_integeral + u_i' * Ru * u_i * simulation_interval;
    x0 = sub_xhost_fixed_interval(end, :)';
end
J_du = norm(diff(u_overall, 1, 1))^2;
J = J + x0'*P_terminal*x0 + J_du;
J = J + B_terminal_cost(x0) - (-0.2); % add terminal cost (SAFETY)
end

%% nonlinear constraint function
% nonlinear constraint function
function [cineq, ceq] = nonlinear_horizon(subsys, u_host, u_adj, x0, x_barrier, d_barrier, N, ORDER)
% suppose ||x0-x_barrier||>=d_barrier
global sample_interval Nc
if ~exist('N', 'var')
    N = 1;
end

if ~exist('ORDER', 'var')
    ORDER = 1;
end

if ORDER == 1
    u_overall = [u_host u_adj];
elseif ORDER == 2
    u_overall = [u_adj(:,1) u_host u_adj(:,2)];
else
    u_overall = [u_adj u_host];
end

cineq = zeros(N, 1); ceq = [];
% U = [u_host u_adj]';
for i = 1:N
    u_i = u_overall(min(i, Nc), :)';
    [~, x0_new] = ode45(@(t, x) subsys(x, u_i), [0 sample_interval], x0);
    x0_new = x0_new(end, :)';
    cineq(i, :) = d_barrier - norm(x0_new([3, 5])-x_barrier);
    x0 = x0_new; % update state
end
end

function [cineq, ceq] = nonlinear_gamma_horizon(subsys, gamma, u_host, u_adj, u_host_pre, u_adj_pre, x0, x_barrier, d_barrier, N, ORDER)
% suppose ||x0-x_barrier||>=d_barrier
global sample_interval Nc
if ~exist('N', 'var')
    N = 1;
end

if ~exist('ORDER', 'var')
    ORDER = 1;
end

cineq = zeros(N, 1); ceq = [];
u_host = gamma(1)*u_host + (1-gamma(1))*u_host_pre;
u_adj = [gamma(2)*u_adj(:, 1) + (1-gamma(2))*u_adj_pre(:, 1) gamma(3)*u_adj(:, 2) + (1-gamma(3))*u_adj_pre(:, 2)];

if ORDER == 1
    u_overall = [u_host u_adj];
elseif ORDER == 2
    u_overall = [u_adj(:,1) u_host u_adj(:,2)];
else
    u_overall = [u_adj u_host];
end

% U = [u_host u_adj]';
for i = 1:N
    u_i = u_overall(min(i, Nc), :)';
    [~, x0_new] = ode45(@(t, x) subsys(x, u_i), [0 sample_interval], x0);
    x0_new = x0_new(end, :)';
    cineq(i, :) = d_barrier - norm(x0_new([3, 5])-x_barrier);
    x0 = x0_new; % update state
end
end