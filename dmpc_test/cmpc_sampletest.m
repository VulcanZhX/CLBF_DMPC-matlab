clear;

% system definition: x' = -x + u
a = -1; b = 1;

global func
global sample_interval simulation_interval
global steps_per_interval

func = @(x,u) a*x + b*u;

% initial state and target state
x0 = 0.9; xs = 1;

% time parameters

sample_interval = 0.1; simulation_interval = 1e-4;
steps_per_interval = round(sample_interval/simulation_interval); % number of steps in one sample interval
sim_steps = 20; % number of simulation steps

% controller parameters
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = 1; % state cost
R = 1; % control cost

% cost function in prediction horizon
x_log = []; u_log = [];
J_all = cost_function([1 2], x0, xs, Nc, Np, Q, R);

for i_sim = 1:sim_steps
    U_initial = [0.2 0.1];
    % define optimization problem
    lb = [-1 -1]; % lower bound
    ub = [2 2]; % upper bound
    opt_problem = createOptimProblem('fmincon','objective',@(U)cost_function(U, x0, xs, Nc, Np, Q, R),...
                'lb',lb,'ub',ub,'x0',U_initial);

    [U_opt, J_opt] = fmincon(opt_problem);

    % next state under U_opt(1)
    [t,x_sol]=ode45(@(t,x) func(x,U_opt(1)), [0 sample_interval], x0);
    x0_new = x_sol(end); % new state
    x_log = [x_log; x0]; % log the state
    u_log = [u_log; U_opt(1)]; % log the control input
    x0 = x0_new; % update state
end

plot(1:sim_steps, x_log, 'r-');
figure;
plot(1:sim_steps, u_log, 'b-');


function J = cost_function(U, x0, xs, Nc, Np, Q, R)
    global func
    global sample_interval simulation_interval
    global steps_per_interval

    J = 0;
    for i = 1:Np
        if i <= Nc
            u_i = U(i);
        else
            u_i = U(Nc);  % Nc <= Np
        end
        [t,x_i]=ode45(@(t,x) func(x,u_i), [0 sample_interval], x0);
        x_i_fixed_interval = interp1(t, x_i, (0:simulation_interval:sample_interval));
        % calculate cost in one sample interval
        discreted_normed_integeral = sum(trace(norm(x_i_fixed_interval-xs)^2))/steps_per_interval;
        J = J + Q*discreted_normed_integeral;
    end
end