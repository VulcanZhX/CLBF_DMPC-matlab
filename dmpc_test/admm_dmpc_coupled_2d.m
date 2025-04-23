%% DMPC Implementation for Coupled 2D Systems via ADMM
clear
% system definition: x' = Ax + Bu (2D system with 2D input)
A = [0.2 -0.41; 0 0.2]; B = diag([0.5 1]);
global func
global Nc Np Q Ru
% global sample_interval simulation_interval
func = @(x,u) A*x  + B*u;
func1 = @(x_local,u,v_adj) A(1, 1)*x_local + B(1, :)*u + v_adj;
func2 = @(x_local,u,v_adj) A(2, 2)*x_local + B(2, :)*u + v_adj;
%func_uv = @(x, u, v) A*x + B*u + v; % coupled system function

% initial state and target state
x0_initial = [-0.2 0.5]'; xs = [0 0]'; x0 = x0_initial;
% time parameters
% sample_interval = 0.1; simulation_interval = 1e-3;
sim_steps = 300; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 10; % number of iterations within one sample interval

% controller parameters
Nc = 2; % controller horizon
Np = 2; % prediction horizon
Q = eye(2); % state cost
Ru = 1e-3*eye(2); % control cost

% log the state and control input
x_log = []; u_log = [];

% initial input guess and bounds
U_initial = zeros(Nc, 2);
U_lb = -1*ones(Nc, 2); % lower bound
U_ub = 1*ones(Nc, 2); % upper bound

%% CODE INSPECTION
v1_initial = [1 1]'; % initial guess for v1(con: v1=A12x2)
V_initial = ones(Nc, 2); % same size (Nc*n_in) as U_initial
UV_init = [U_initial(:, 1) v1_initial]; % initial guess for U_host and v (combined for opt.)
localcost_s1 = @(UV)cost_function(func1, UV, U_initial(:, 2), x0(1), xs(1), 1);

opt_options = optimoptions('fmincon', 'Display', 'iter');
opt_problem1 = createOptimProblem('fmincon', 'objective', localcost_s1, ...
    'x0', UV_init, 'lb', U_lb, 'ub', U_ub, 'options', opt_options);
[UV_opt1, fval1, exitflag1, info1] = fmincon(opt_problem1);
U_opt1 = UV_opt1(:, 1);
V_opt1 = UV_opt1(:, 2);
dual_all = zeros(Np, 2);

%% It's really exhausting to figure out the correct shape of each matrix and dual variables in the Lagrangian Function!
% [UV_init(:, 2)' repmat(UV_init(2, 2), 1, Np - Nc)] - A(1, 2)*forward_simulation(func1, UV_init, U_initial(:, 2), x0(1), 1);
lagrange_localcost_s1 = @(UV, dual_all) localcost_s1(UV) + trace(dual_all(:, 1)*[UV_init(:, 2)' repmat(UV_init(2, 2), 1, Np - Nc)]) ...
    - trace(dual_all(:, 2)*A(1, 2)*forward_simulation(func1, UV, U_initial(:, 2), x0(1), 1));
lagrange_localcost_s1(UV_init, dual_all)

opt_options = optimoptions('fmincon', 'Display', 'off');
dual_overall = ones(Np, 2); % initial guess for dual variables(row: time step, column: subsystem label)

%% Rolling Optimization
for i_sim = 1:sim_steps
    tic
    % optimization problem for s1
    U_host_init = U_initial(:, 1);
    U_adj_init = U_initial(:, 2);
    UV_init = [U_host_init, V_initial(:, 1)];
    % local cost function
    localcost_s1 = @(UV)cost_function(func1, UV, U_adj_init, x0(1), xs(1), 1);
    % local Lagrangian cost function
    lagrange_localcost_s1 = @(UV, dual_all) localcost_s1(UV) + trace(dual_all(:, 1)*[UV(:, 2)' repmat(UV(2, 2), 1, Np - Nc)]) ...
                            - trace(dual_all(:, 2)*A(2, 1)*forward_simulation(func1, UV, U_initial(:, 2), x0(1), 1));
    % optimizing primal variables of s1 Lagrangian cost function
    primal_opt_problem1 = createOptimProblem('fmincon', 'objective', @(UV)lagrange_localcost_s1(UV, dual_overall), ...
        'x0', UV_init, 'lb', U_lb, 'ub', U_ub, 'options', opt_options);
    [UV_opt1, lagrange_cost1] = fmincon(primal_opt_problem1);

    % optimization problem for s2
    U_host_init = U_initial(:, 2);
    U_adj_init = U_initial(:, 1);
    UV_init = [U_host_init, V_initial(:, 2)];
    % local cost function
    localcost_s2 = @(UV)cost_function(func2, UV, U_adj_init, x0(2), xs(2), 2);
    % local Lagrangian cost function
    lagrange_localcost_s2 = @(UV, dual_all) localcost_s2(UV) + trace(dual_all(:, 2)*[UV(:, 2)' repmat(UV(2, 2), 1, Np - Nc)]) ...
        - trace(dual_all(:, 1)*A(1, 2)*forward_simulation(func2, UV, U_initial(:, 1), x0(2), 2));
    % optimizing primal variables of s2 Lagrangian cost function
    primal_opt_problem2 = createOptimProblem('fmincon', 'objective', @(UV)lagrange_localcost_s2(UV, dual_overall), ...
        'x0', UV_init, 'lb', U_lb, 'ub', U_ub, 'options', opt_options);
    [UV_opt2, lagrange_cost2] = fmincon(primal_opt_problem2);

    % optimizing dual variables
    % dual variables for s1
    % dual_opt_option = optimoptions('fmincon', 'Display', 'notify', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true);
    % step size for dual variable update
    step1 = 0.01; 
    step2 = 0.01;
    % dual variable update
    UV_all = [UV_opt1(:, 1) UV_opt2(:, 1) UV_opt1(:, 2) UV_opt2(:, 2)]; % [U1 U2 V1 V2]
    dual_overall(:, 1) = dual_overall(:, 1) + step1*(UV_all(:, 3) - A(1, 2)*(forward_simulation(func2, UV_opt2, U_initial(:, 1), x0(2), 2))');
    dual_overall(:, 2) = dual_overall(:, 2) + step2*(UV_all(:, 4) - A(2, 1)*(forward_simulation(func1, UV_opt1, U_initial(:, 2), x0(1), 1))');

    % update
    U_initial = UV_all(:, 1:2);
    V_initial = UV_all(:, 3:4);
    x_new = func(x0, U_initial(:, 1));
    x_log = [x_log x0]; % log state
    u_log = [u_log U_initial(:, 1)']; % log control input
    x0 = x_new; % update state

    % display
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f\n', ...
        i_sim, toc, norm(x0 - xs));
end

%% forward simulation
function x_sample = forward_simulation(func, u_host_v, u_adj, x0_local, ORDER)
global Np Nc
x_sample = x0_local; % initialize state
u_host = u_host_v(:, 1);
v = u_host_v(:, 2); % extract host input and adjoint state

if ~exist('ORDER', 'var')
    ORDER = 1; % default order
end

if ORDER == 1
    u_overall = [u_host u_adj];
else
    u_overall = [u_adj u_host];
end

for i = 1:Np-1
    u_i = u_overall(min(i, Nc), :)';
    v_i = v(min(i, Nc), :)';
    x0_local_new = func(x0_local, u_i, v_i);
    x_sample = [x_sample x0_local_new]; % update state
    x0_local = x0_local_new;
end
end

%% local cost function with external v
function J_local = cost_function(subsys, u_host_v, u_adj, x0_local, xs_local, ORDER)
% calculate cost of subsystem S_i
global Nc Np Q Ru
u_host = u_host_v(:, 1);
v = u_host_v(:, 2); % extract host input and adjoint state

J_local = 0; % initialize cost
if ~exist('ORDER', 'var')
    ORDER = 1; % default order
end
if ORDER == 1
    u_overall = [u_host u_adj];
else
    u_overall = [u_adj u_host];
end

for i = 1:Np
    u_i = u_overall(min(i, Nc), :)';
    v_i = v(min(i, Nc), :)';
    J_local = J_local + (x0_local-xs_local)' * Q(ORDER, ORDER) * (x0_local-xs_local) + u_i(ORDER)' * Ru(ORDER, ORDER) * u_i(ORDER); % cost update
    x0_new_local = subsys(x0_local, u_i, v_i);
    x0_local = x0_new_local; % update state
end
end
