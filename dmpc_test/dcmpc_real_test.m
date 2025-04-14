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

loc_func1 = @(x,u)[
    (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/x(3))*x(1);
    (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/x(3))*x(1)-kB*exp(-EBR/x(3))*x(2);
    (Ff1*T0+FR*x(9)-F1*x(3))/V1-(kA*exp(-EAR/x(3))*x(1)*dHA+kB*exp(-EBR/x(3))*x(2)*dHB)/Cp+u(1)/(Cp*V1);
    ];

loc_func2 = @(x,u)[
    (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/x(6))*x(4);
    (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/x(6))*x(4)-kB*exp(-EBR/x(6))*x(5);
    (Ff2*T0+F1*x(3)-F2*x(6))/V2-(kA*exp(-EAR/x(6))*x(4)*dHA+kB*exp(-EBR/x(6))*x(5)*dHB)/Cp+u(2)/(Cp*V2);
    ];

loc_func3 = @(x,u)[
    (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    (F2*x(6)-(FD+FR)*x(9)-F3*x(9))/V3+u(3)/(Cp*V3);
    ];

% initial state and target state
xs1=0;
xs2=0;
xs3=499.479234575489;
xs4=0;
xs5=0;
xs6=475.482353426126;
xs7=0;
xs8=0;
xs9=314.757389283257;
xs = [xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9]';
x0_initial = [0.5 0.5 470 0.2 0.7 450 0.1 0.8 341]'; x0 = x0_initial;

% time parameters
global sample_interval simulation_interval
sample_interval = 0.1; simulation_interval = 1e-4;
sim_steps = 20000; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 1; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 4; % controller horizon
Np = 6; % prediction horizon
Q = diag([1, 1, 2, 1, 1, 2, 1, 1, 2]); % state cost
Ru = 1e-4*eye(3); % control cost

% initial input guess and bounds
U0 = 12.6e5 * ones(Nc, 3); % initial guess
U_lb = (12.6e5 - 5e5)*ones(Nc, 3); % lower bound
U_ub = (12.6e5 + 5e5)*ones(Nc, 3); % upper bound

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
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 1);
        opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'options', opt_options);
        [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input 2
        U_host_init = U0(:, 2);
        U_adj_init = [U0(:, 1) U0(:, 3)];
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 2);
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 2), 'ub', U_ub(:, 2), 'options', opt_options);
        [U_opt2, J_opt2] = fmincon(opt_problem2);

        % optimization problem for input 3
        U_host_init = U0(:, 3);
        U_adj_init = [U0(:, 1) U0(:, 2)];
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 3);
        opt_problem3 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 3), 'ub', U_ub(:, 3), 'options', opt_options);
        [U_opt3, J_opt3] = fmincon(opt_problem3);

        % convex combination of the three optimization problems
        weightcost_obj = @(gamma) sub_cost_function(func, gamma(1)*U_opt1 + (1-gamma(1))*U0(:, 1),...
            [gamma(2)*U_opt2 + (1-gamma(2))*U0(:, 2), gamma(3)*U_opt3 + (1-gamma(3))*U0(:, 3)], x0, xs, 1); % gamma in [0, 1]

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
    if norm(x0 - xs) < 25
        sample_interval = 0.005;
        % R = 1e-2*eye(3);
    end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0 - xs) < 10
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end


% function x0_new = sysfwd(sys, x0, u)
% % sysfwd: forward simulation of the system
% global sample_interval
% [~, x_sol] = ode45(@(t,x) sys(x, u), [0 sample_interval], x0);
% x0_new = x_sol(end, :)';
% end

% function J = sub_cost_function(subsys, u_host, u_adj, x0, xs)
% % sub_cost_function: cost function for the subsystem
% global sample_interval simulation_interval Nc Np Q Ru

% % remember to delete this
% x0_host_s = xs; % target state

% J = 0; % initialize cost

% u_overall = [u_host u_adj];

% for i = 1:Np
%     if i <= Nc
%         u_i = u_overall(i, :)';
%     else
%         u_i = u_overall(Nc, :)';  % Nc <= Np
%     end
%     [t, sub_xhost_sol] = ode45(@(t, x) subsys(x, u_i), [0 sample_interval], x0);
%     sub_xhost_fixed_interval = interp1(t, sub_xhost_sol, (0:simulation_interval:sample_interval-simulation_interval));

%     % calculate cost in one sample interval
%     delta_subx_array = sub_xhost_fixed_interval - repmat(x0_host_s', size(sub_xhost_fixed_interval, 1), 1); % X - Xs
%     sub_weighted_norm = trace(delta_subx_array * Q * delta_subx_array');
%     sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
%     J = J + sub_discreted_normed_integeral + u_i' * Ru * u_i * simulation_interval;
%     x0 = sub_xhost_fixed_interval(end, :)';
% end
% end

function x_nxt = sysfwd(sys, x0, u)
% sysfwd: forward simulation of the system (one sample interval)
global sample_interval simulation_interval
x_nxt = zeros(size(x0, 1), round(sample_interval/simulation_interval));
fwd_step = round(sample_interval/simulation_interval);
for i = 1:fwd_step
    eqn = sys(x0, u);
    x0 = x0 + eqn*simulation_interval; % update state
    x_nxt(:, i) = x0; % store state
end
end

function J = sub_cost_function(subsys, u_host, u_adj, x0, xs, ORDER)
global sample_interval simulation_interval Nc Np Q Ru
J = 0; % initialize cost
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
    x_sample = sysfwd(subsys, x0, u_i); % size: size(x0, 2) x round(sample_interval/simulation_interval)
    % calculate cost in one sample interval
    delta_subx_array = x_sample - repmat(xs, 1, size(x_sample, 2)); % X - Xs
    sub_weighted_norm = trace(delta_subx_array' * Q * delta_subx_array);
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i' * Ru * u_i * sample_interval;
    x0 = x_sample(:, end); % update state
end
end