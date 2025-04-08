clear;

% system definition: x' = Ax + Bu (2D system with 2D input)
A = -diag([1 2]); B = diag([0.8 0.5]);

% global func1 func2
global sample_interval simulation_interval
global steps_per_interval

func1 = @(x1, x2, u1, u2) A(1,:)*[x1; x2]  + B(1,:)*[u1; u2];
func2 = @(x,u) A(2,:)*[x1; x2]  + B(2,:)*[u1; u2];

% initial state and target state
x0_initial = [1 1]'; xs = [0.5 0]'; x0 = x0_initial;

% time parameters
sample_interval = 0.1; simulation_interval = 1e-4;
steps_per_interval = round(sample_interval/simulation_interval); % number of steps in one sample interval
sim_steps = 10; % number of simulation steps, overall time = sim_steps*sample_interval

% controller parameters
Nc = 2; % controller horizon
Np = 5; % prediction horizon
Q = eye(2); % state cost
R = 1e-4*eye(2); % control cost

% cost function in prediction horizon
x_log = []; u_log = [];
u_initial = [0.1 0.1; 0.3 0.2];

sub_x0_struct = struct('xhost', x0(1), 'xhost_s', xs(1), 'xadj', x0(2), 'xadj_s', xs(2));
sub_U_struct = struct('uhost', u_initial(:, 1), 'uadj', u_initial(:, 2));

J_sys1 = sub_cost_function(func1, sub_U_struct, sub_x0_struct, Nc, Np, Q(1, 1), R(1, 1));
disp(['Cost for system 1: ', num2str(J_sys1)]);


function J = sub_cost_function(subsys, sub_U_struct, sub_x0_struct, Nc, Np, sub_Q, sub_R)
global sample_interval simulation_interval
global steps_per_interval

J = 0;

x0_host = sub_x0_struct.xhost;
x0_host_s = sub_x0_struct.xhost_s;
x0_adj = sub_x0_struct.xadj;
u_adj = sub_U_struct.uadj;

for i = 1:Np
    if i <= Nc
        u_i = sub_U_struct.uhost(i, :)';
        u_adj = sub_U_struct.uadj(i, :)';
    else
        u_i = sub_U_struct.uhost(Nc, :)';  % Nc <= Np
        u_adj = sub_U_struct.uadj(Nc, :)';
    end
    [t, sub_xhost_sol] = ode45(@(t, x) subsys(x, x0_adj, u_i, u_adj), [0 sample_interval], x0_host);
    sub_xhost_fixed_interval = interp1(t, sub_xhost_sol, (0:simulation_interval:sample_interval-simulation_interval))';

    % calculate cost in one sample interval
    delta_subx_array = sub_xhost_fixed_interval - repmat(x0_host_s', size(sub_xhost_fixed_interval, 1), 1); % X - Xs
    sub_weighted_norm = trace(delta_subx_array' * sub_Q * delta_subx_array);
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i' * sub_R * u_i * simulation_interval;
    x0_host = sub_xhost_fixed_interval(end, :);
end
end