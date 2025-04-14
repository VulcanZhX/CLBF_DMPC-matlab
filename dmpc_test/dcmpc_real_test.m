clear
rng('default')

% system parameters
global F1 F2 F3 FD FR Ff1 Ff2 V1 V2 V3 alA alB alC kA kB EAR EBR dHA dHB Cp T0 xA0
F1=35.5; F2=43.5; F3=15.5;
FD=0.504; FR=50.4;
Ff1=5; Ff2=5;
V1=1*1000; V2=0.5*1000; V3=0.012*1000;
alA=3.5; alB=1;alC=0.5;
kA=2.77e3*3600; kB=2.5e3*3600;
EA=50000; EB=60000;
R=8.314; EAR=EA/R; EBR=EB/R;
MW=250e-3; dHA=-60000/MW; dHB=-70000/MW;
Cp=4.2e3; T0=313; xA0=1;

% system definition
func = @(x,u) [
    (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/x(3))*x(1);
    (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/x(3))*x(1)-kB*exp(-EBR/x(3))*x(2);
    (Ff1*T0+FR*x(9)-F1*x(3))/V1-(kA*exp(-EAR/x(3))*x(1)*dHA+kB*exp(-EBR/x(3))*x(2)*dHB)/Cp+u(1)/(Cp*V1);
    (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/x(6))*x(4);
    (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/x(6))*x(4)-kB*exp(-EBR/x(6))*x(5);
    (Ff2*T0+F1*x(3)-F2*x(6))/V2-(kA*exp(-EAR/x(6))*x(4)*dHA+kB*exp(-EBR/x(6))*x(5)*dHB)/Cp+u(2)/(Cp*V2);
    (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    (F2*x(6)-(FD+FR)*x(9)-F3*x(9))/V3+u(3)/(Cp*V3);
    ];

% initial state and target state
x0_initial = [0.2 0.5 496 0.4 0.5 489 0.25 0.6 350]'; xs = [0.383 0.581 450 0.391 0.572 445 0.172 0.748 449]'; x0 = x0_initial;

% time parameters
global sample_interval simulation_interval
sample_interval = 0.05; simulation_interval = 1e-4;
sim_steps = 2000; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 50; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = diag([1, 1, 1, 1, 1, 1, 1, 1, 1]); % state cost
Ru = 1e-8*eye(3); % control cost

% initial input guess and bounds
U0 = 12.6e5 * ones(Nc, 3); % initial guess
U_lb = (12.6e5 - 4e5)*ones(Nc, 3); % lower bound
U_ub = (12.6e5 +4e5)*ones(Nc, 3); % upper bound

% log the state and control input
x_log = []; u_log = [];

% rolling optimization
opt_options = optimoptions('fmincon', 'Display', 'off');
simulation_timer = 0;
for i_sim = 1:sim_steps
    tic; % begin current loop timer
    for i_p = 1:max_p_iteration
        % optimization problem for input 1
        U_host_init = U0(:, 1);
        U_adj_init = U0(:, 2:3);
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs);
        opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'options', opt_options);
        [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input 2
        U_host_init = U0(:, 2);
        U_adj_init = [U0(:, 1) U0(:, 3)];
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 2), 'ub', U_ub(:, 2), 'options', opt_options);
        [U_opt2, J_opt2] = fmincon(opt_problem2);

        % optimization problem for input 3
        U_host_init = U0(:, 3);
        U_adj_init = [U0(:, 1) U0(:, 2)];
        opt_problem3 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 3), 'ub', U_ub(:, 3), 'options', opt_options);
        [U_opt3, J_opt3] = fmincon(opt_problem3);

        % convex combination of the three optimization problems
        weightcost_obj = @(gamma) sub_cost_function(func, gamma(1)*U_opt1 + (1-gamma(1))*U0(:, 1),...
            [gamma(2)*U_opt2 + (1-gamma(2))*U0(:, 2), gamma(3)*U_opt3 + (1-gamma(3))*U0(:, 3)], x0, xs); % gamma in [0, 1]

        % optimization problem for the convex combination
        opt_problem_comb = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
            'x0', [0.5 0.5 0.5], 'lb', [0 0 0], 'ub', [1 1 1], 'options', opt_options);
        [gamma_opt, J_opt_comb] = fmincon(opt_problem_comb);
        % update control input
        U_opt1_new = gamma_opt(1)*U_opt1 + (1-gamma_opt(1))*U0(:, 1);
        U_opt2_new = gamma_opt(2)*U_opt2 + (1-gamma_opt(2))*U0(:, 2);
        U_opt3_new = gamma_opt(3)*U_opt3 + (1-gamma_opt(3))*U0(:, 3);

        % iteration terminal condition
        if norm([U_opt1_new; U_opt2_new; U_opt3_new] - [U0(:, 1); U0(:, 2); U0(:, 3)]) < 1e-2
            break;
        end
        U0 = [U_opt1_new U_opt2_new U_opt3_new]; % update initial guess
    end
    % update state
    [t, x_global_sol] = ode45(@(t, x) func(x, U0(1, :)'), [0 sample_interval], x0);
    x0_new = x_global_sol(end, :)'; % new state
    x_log = [x_log x0_new]; % log the state
    u_log = [u_log U0(1, :)']; % log the control input
    x0 = x0_new; % update state
    % display elapsed time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f\n', i_sim, toc, norm(x0 - xs));
    % terminal condition
    if norm(x0 - xs) < 15
        sample_interval = 0.005;
        % R = 1e-2*eye(3);
    end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0 - xs) < 1
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end


function x0_new = sysfwd(sys, x0, u)
% sysfwd: forward simulation of the system
global sample_interval
[~, x_sol] = ode45(@(t,x) sys(x, u), [0 sample_interval], x0);
x0_new = x_sol(end, :)';
end

function J = sub_cost_function(subsys, u_host, u_adj, x0, xs)
% sub_cost_function: cost function for the subsystem
global sample_interval simulation_interval Nc Np Q Ru

% remember to delete this
x0_host_s = xs; % target state

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
    J = J + sub_discreted_normed_integeral + u_i' * Ru * u_i * simulation_interval;
    x0 = sub_xhost_fixed_interval(end, :)';
end
end