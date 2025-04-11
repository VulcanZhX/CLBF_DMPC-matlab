clear;

% system definition: x' = Ax + Bu (2D system with 2D input)
A = [-1 0.3; -0.15 -1.5]; B = [0.6 0.1;0.2 0.8];

% global func1 func2
global sample_interval simulation_interval
global steps_per_interval

func = @(x, u) A*x + B*u; % system function
func1 = @(x1, x2, u1, u2) A(1,:)*[x1; x2]  + B(1,:)*[u1; u2];
func2 = @(x1, x2, u1, u2) A(2,:)*[x1; x2]  + B(2,:)*[u1; u2];

% initial state and target state
x0_initial = [1.5 3]'; xs = [0 0]'; x0 = x0_initial;

% time parameters
sample_interval = 0.05; simulation_interval = 1e-4;
steps_per_interval = round(sample_interval/simulation_interval); % number of steps in one sample interval
sim_steps = 100; % number of simulation steps, overall time = sim_steps*sample_interval

% controller parameters
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = eye(2); % state cost
R = 1e-3*eye(2); % control cost

% cost function in prediction horizon
x_log = []; u_log = [];
u_initial = [0.1 0.1; 0.3 0.2];

sub_x0_struct = struct('xhost', x0(1), 'xhost_s', xs(1), 'xadj', x0(2), 'xadj_s', xs(2));
sub_U_struct = struct('uhost', u_initial(:, 1), 'uadj', u_initial(:, 2));

J_sys1 = sub_cost_function(func1, [u_initial(:, 1); 1000; 1000], u_initial(:, 2), sub_x0_struct, Nc, Np, Q(1, 1), R(1, 1));

U_initial = [0.2 0.1; 0.1 0];
lb = [-3 -3; -3 -3];
ub = [3 3; 3 3];
options = optimoptions('fmincon', 'Display', 'off');
tic;
for i_sim = 1:sim_steps
    % optimization problem for system 1
    sub_x0_struct = struct('xhost', x0(1), 'xhost_s', xs(1), 'xadj', x0(2), 'xadj_s', xs(2));
    opt_problem_1 = createOptimProblem('fmincon', 'objective', @(u) sub_cost_function(func1, u, u_initial(:, 2), sub_x0_struct, Nc, Np, Q(1, 1), R(1, 1)), ...
        'x0', U_initial(:), 'lb', lb(:), 'ub', ub(:), 'options', options);
    [U_opt1, J_opt1] = fmincon(opt_problem_1); 
    % optimization problem for system 2
    sub_x0_struct = struct('xhost', x0(2), 'xhost_s', xs(2), 'xadj', x0(1), 'xadj_s', xs(1));
    opt_problem_2 = createOptimProblem('fmincon', 'objective', @(u) sub_cost_function(func2, u, u_initial(: ,1), sub_x0_struct, Nc, Np, Q(2, 2), R(2, 2)), ...
        'x0', U_initial(:), 'lb', lb(:), 'ub', ub(:), 'options', options);
    [U_opt2, J_opt2] = fmincon(opt_problem_2);

    % update control input and next state
    [t, x_full_sol] = ode45(@(t, x) func(x, [U_opt1(1,:), U_opt2(1,:)]'), [0 sample_interval], x0);
    x0_new = x_full_sol(end, :)';
    x_log = [x_log; x0_new'];
    u_log = [u_log; [U_opt1(1, :), U_opt2(1, :)]];
    x0 = x0_new; 
    u_initial = [U_opt1 U_opt2];
    u_initial = u_initial(1:Nc,:);
end
plot(x_log(:,1), x_log(:,2), '-*')
toc;

function J = sub_cost_function(subsys, u_host, u_adj, sub_x0_struct, Nc, Np, sub_Q, sub_R)
global sample_interval simulation_interval
global steps_per_interval

J = 0;

x0_host = sub_x0_struct.xhost;
x0_host_s = sub_x0_struct.xhost_s;
x0_adj = sub_x0_struct.xadj;
% u_adj = sub_U_struct.uadj;

for i = 1:Np
    if i <= Nc
        u_i = u_host(i, :)';
        u_j = u_adj(i, :)';
    else
        u_i = u_host(Nc, :)';  % Nc <= Np
        u_j = u_adj(Nc, :)';
    end
    [t, sub_xhost_sol] = ode45(@(t, x) subsys(x, x0_adj, u_i, u_j), [0 sample_interval], x0_host);
    sub_xhost_fixed_interval = interp1(t, sub_xhost_sol, (0:simulation_interval:sample_interval-simulation_interval))';

    % calculate cost in one sample interval
    delta_subx_array = sub_xhost_fixed_interval - repmat(x0_host_s', size(sub_xhost_fixed_interval, 1), 1); % X - Xs
    sub_weighted_norm = trace(delta_subx_array' * sub_Q * delta_subx_array);
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i' * sub_R * u_i * simulation_interval;
    x0_host = sub_xhost_fixed_interval(end, :);
end
end