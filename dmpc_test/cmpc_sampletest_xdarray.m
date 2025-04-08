clear;

% system definition: x' = Ax + Bu (2D system with 2D input)
A = -diag([1 2]); B = diag([0.8 0.5]);

global func
global sample_interval simulation_interval
global steps_per_interval

func = @(x,u) A*x  + B*u;

% initial state and target state
x0_initial = [1 1]'; xs = [0.5 0]'; x0 = x0_initial;

% time parameters
sample_interval = 0.1; simulation_interval = 1e-4;
steps_per_interval = round(sample_interval/simulation_interval); % number of steps in one sample interval
sim_steps = 10; % number of simulation steps, overall time = sim_steps*sample_interval

% controller parameters
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = eye(2); % state cost
R = 1e-4*eye(2); % control cost

% cost function in prediction horizon
x_log = []; u_log = [];
J_all = cost_function([1 1; 1 1], x0, xs, Nc, Np, Q, R);

for i_sim = 1:sim_steps
    U_initial = [0.2 0.1; 0.1 0];
    % define optimization problem
    lb = [-2 -2; -2 -2]; % lower bound
    ub = [4 4; 4 4]; % upper bound
    opt_problem = createOptimProblem('fmincon','objective',@(U)cost_function(U, x0, xs, Nc, Np, Q, R),...
        'lb',lb,'ub',ub,'x0',U_initial);

    [U_opt, J_opt] = fmincon(opt_problem);

    % next state under U_opt(1)
    [t,x_sol]=ode45(@(t,x) func(x,U_opt(1, :)'), [0 sample_interval], x0);
    x0_new = x_sol(end,:)'; % new state
    x_log = [x_log x0]; % log the state
    u_log = [u_log U_opt(1,:)']; % log the control input
    x0 = x0_new; % update state
end



syms t;
uf(t) = [0; 0]*t;
for i = 1:length(u_log)
    uf(t) = uf(t) + u_log(:, i)*(heaviside(t-(i-1)*sample_interval) - heaviside(t-i*sample_interval));
end

uf_handle = matlabFunction(uf);
% [t, x_sol]=ode45(@(t, x)func(x, uf_handle(t)), [0 sim_steps*sample_interval], x0_initial);
% figure
% subplot(2,1,1);
% plot(t, x_sol);
% title('State evolution');
% xlabel('Time (s)');
% ylabel('State x');
% subplot(2,1,2);
% t = 0:simulation_interval:sim_steps*sample_interval;
% plot(t, uf_handle(t));
% title('Control input evolution');
% xlabel('Time (s)');
% ylabel('Control input u');




function J = cost_function(U, x0, xs, Nc, Np, Q, R)
global func
global sample_interval simulation_interval
global steps_per_interval

J = 0;
for i = 1:Np
    if i <= Nc
        u_i = U(i, :)';
    else
        u_i = U(Nc, :)';  % Nc <= Np
    end
    [t,x_i]=ode45(@(t,x) func(x,u_i), [0 sample_interval], x0);
    x_i_fixed_interval = interp1(t, x_i, (0:simulation_interval:sample_interval-simulation_interval));

    % calculate cost in one sample interval
    delta_x_array = x_i_fixed_interval-repmat(xs', size(x_i_fixed_interval, 1), 1); %X - Xs
    weighted_norm = trace(delta_x_array*Q*delta_x_array');
    discreted_normed_integeral = weighted_norm*simulation_interval;
    J = J + discreted_normed_integeral + u_i'*R*u_i*sample_interval;
    x0 = x_i_fixed_interval(end,:);
end
end