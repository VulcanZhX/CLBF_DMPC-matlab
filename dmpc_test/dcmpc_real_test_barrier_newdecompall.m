clear
rng(12345)

% system parameters
global F1 F2 F3 FD FR Ff1 Ff2 V1 V2 V3 alA alB alC kA kB EAR EBR dHA dHB Cp T0 xA0
F1=74.53; F2=75.03; F3=8.168;
FD=0.662; FR=66.2;
Ff1=8.33; Ff2=0.5;
V1=13.41; V2=13.5; V3=0.5;
alA=3.5; alB=1.1;alC=0.5;
kA=0.2; kB=0.018;
EAR=-100; EBR=-150;
dHA=-40; dHB=-50;
Cp=2.5; T0=313; xA0=1;

% system definition
func = @(x,u) [
    (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/(x(3)*1))*x(1);
    (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/(x(3)*1))*x(1)-kB*exp(-EBR/(x(3)*1))*x(2);
    ((Ff1*T0+FR*(x(9)*1)-F1*(x(3)*1))/V1-(kA*exp(-EAR/(x(3)*1))*x(1)*dHA+kB*exp(-EBR/(x(3)*1))*x(2)*dHB)/Cp+u(1)/(Cp*V1));
    (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/(x(6)*1))*x(4);
    (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/(x(6)*1))*x(4)-kB*exp(-EBR/(x(6)*1))*x(5);
    (Ff2*T0+F1*(x(3)*1)-F2*(x(6)*1))/V2-(kA*exp(-EAR/(x(6)*1))*x(4)*dHA+kB*exp(-EBR/(x(6)*1))*x(5)*dHB)/Cp+u(2)/(Cp*V2);
    (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    (F2*(x(6)*1)-(FD+FR)*(x(9)*1)-F3*(x(9)*1))/V3+u(3)/(Cp*V3);
    ];

% solve system equilibrium
u_stable = [0 0 0]';
x0_initial = [0.13 0.15 440 0.25 0.1 421 0.12 0.170 404]'; x0 = x0_initial;
[xs, fval] = fsolve(@(x) func(x, u_stable), x0_initial, optimset('Display', 'off'));
% % initial state and target state (backup)
% % xsi solved by fsolve(@(x)func(x,0),x_init)
% xs1=0.696;
% xs2=0.295;
% xs3=322.10;
% xs4=0.66;
% xs5=0.32;
% xs6=322.59;
% xs7=0.41;
% xs8=0.55;
% xs9=322.599;
% xs = [xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9]';
% x0_initial = [0.3 0.5 475 0.2 0.421 450 0.2 0.570 315]'; x0 = x0_initial;

% time parameters
global sample_interval simulation_interval
sample_interval = 0.01; simulation_interval = 1e-4;
sim_steps = 1000; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 10; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 2; % controller horizon
Np = 4; % prediction horizon
Q = diag([1, 1, 0.1, 1, 1, 0.1, 1.2, 1.2, 0.2]); % state cost
Ru = 1e-4*eye(3); % control cost

% initial input guess and bounds
U0 =  zeros(Nc, 3); % initial guess
U_lb = repmat([-1000 -1000 -750], Nc, 1); % lower bound
U_ub = repmat([1000 1000 750], Nc, 1); % upper bound

% log the state and control input
x_log = []; u_log = [];

% barriers
x_barrier = [0.4 0.11]'; % barrier point (for reactor 2)
d_barrier = 0.077; % distance to barrier
global B_terminal_cost
B_terminal_cost = @(x) (d_barrier - norm(x([4, 5])-x_barrier));



%% rolling optimization
opt_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','interior-point','EnableFeasibilityMode',true);
simulation_timer = 0;
for i_sim = 1:sim_steps
    tic; % begin current loop timer
    for i_p = 1:max_p_iteration
        % optimization problem for input 1
        % U_host_init = U0(:, 1);
        % U_adj_init = U0(:, 2:3);
        % localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 1);
        % nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, Np, 1);
        % opt_problem1 = createOptimProblem('fmincon', 'objective', localcost, ...
        %     'x0', U_host_init, 'lb', U_lb(:, 1), 'ub', U_ub(:, 1), 'nonlcon', nonlin_handle, 'options', opt_options);
        % [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input 1,2
        U_host_init = [U0(:, 1) U0(:, 2) U0(:, 3)];
        U_adj_init = [];
        localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 2);
        nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, Np, 2);
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost, ...
            'x0', U_host_init, 'lb', U_lb(:, 1:3), 'ub', U_ub(:, 1:3), 'nonlcon', nonlin_handle, 'options', opt_options);
        [U_opt, J_opt2] = fmincon(opt_problem2);
        U_opt1 = U_opt(:, 1);
        U_opt2 = U_opt(:, 2);
        U_opt3 = U_opt(:, 3);

        % optimization problem for input 3
        % U_host_init = U0(:, 3);
        % U_adj_init = [U0(:, 1) U0(:, 2)];
        % localcost = @(U)sub_cost_function(func, U, U_adj_init, x0, xs, 3);
        % nonlin_handle = @(U)nonlinear_horizon(func, U, U_adj_init, x0, Np, 3);
        % opt_problem3 = createOptimProblem('fmincon', 'objective', localcost, ...
        %     'x0', U_host_init, 'lb', U_lb(:, 3), 'ub', U_ub(:, 3), 'nonlcon', nonlin_handle, 'options', opt_options);
        % [U_opt3, J_opt3] = fmincon(opt_problem3);

        % convex combination of the three optimization problems
        weightcost_obj = @(gamma) sub_cost_function(func, gamma(1)*U_opt1 + (1-gamma(1))*U0(:, 1),...
            [gamma(2)*U_opt2 + (1-gamma(2))*U0(:, 2), gamma(3)*U_opt3 + (1-gamma(3))*U0(:, 3)], x0, xs, 1); % gamma in [0, 1]
        % nonlinear constraint for the convex combination
        nonlin_handle_comb = @(gamma)nonlinear_gamma_horizon(func, gamma, U_opt1, [U_opt2 U_opt3], U0(:, 1), U0(:, 2:3), x0, Np, 1);
        % optimization problem for the convex combination
        opt_problem_comb = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
            'x0', [0.5 0.5 0.5], 'lb', [0 0 0], 'ub', [1 1 1], 'nonlcon', nonlin_handle_comb, 'options', opt_options);
        [gamma_opt, J_opt_comb] = fmincon(opt_problem_comb);
        % update control input
        U_opt1_new = gamma_opt(1)*U_opt1 + (1-gamma_opt(1))*U0(:, 1);
        U_opt2_new = gamma_opt(2)*U_opt2 + (1-gamma_opt(2))*U0(:, 2);
        U_opt3_new = gamma_opt(3)*U_opt3 + (1-gamma_opt(3))*U0(:, 3);
        % update initial guess
        U0 = [U_opt1_new U_opt2_new U_opt3_new]; % update initial guess
        % iteration terminal condition
        if norm([U_opt1_new; U_opt2_new; U_opt3_new] - [U0(:, 1); U0(:, 2); U0(:, 3)]) < 3
            break;
        end
    end
    % update state
    [t, x_global_sol] = ode45(@(t, x) func(x, U0(1, :)'), [0 sample_interval], x0);
    x0_new = x_global_sol(end, :)'; % new state
    x_log = [x_log x0_new]; % log the state
    u_log = [u_log U0(1, :)']; % log the control input
    x0 = x0_new; % update state
    % display elapsed time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f, safety: %.4f, J_cmb: %.4f\n', ...
        i_sim, toc, norm(x0 - xs), B_terminal_cost(x0), J_opt_comb);
    % terminal condition
    % if norm(x0 - xs) < 5
    %     % sample_interval = 0.01;
    %     R = 5e-3*eye(3);
    % end
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0([1:2, 4:5, 7:8]) - xs([1:2, 4:5, 7:8])) < 0.04 && norm(x0([3, 6, 9])-xs([3, 6, 9])) < 1
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
    % if norm(x0(4:5) - x_barrier) < 0.15
end

x_log = [x0_initial x_log]; % add initial state to log

% simulation time
fprintf('Total simulation time: %.4f seconds\n', simulation_timer);
% plot temperature
figure
grid on
xlabel('sim_step')
plot(x_log(3,:),'linewidth',1.5)
hold on
plot(x_log(6,:),'linewidth',1.5)
hold on
plot(x_log(9,:),'linewidth',1.5)
ylabel('Temperature (K)')
legend('T1', 'T2', 'T3')

% plot other states
figure
grid on
plot(x_log(1,:), x_log(2,:), 'linewidth', 1, 'Color', 'r')
hold on
plot(x_log(4,:), x_log(5,:), 'linewidth', 1, 'Color', 'g')
hold on
plot(x_log(7,:), x_log(8,:), 'linewidth', 1, 'Color', 'b')
hold on
scatter(xs(1), xs(2), 'r*')
hold on
scatter(xs(4), xs(5), 'g*')
hold on
scatter(xs(7), xs(8), 'b*')
hold on
fimplicit(@(x1,x2) (x1-x_barrier(1))^2+(x2-x_barrier(2))^2-0.95*d_barrier^2, "--")
xlabel('Concentration of A')
ylabel('Concentration of B')
legend('Reactor 1', 'Reactor 2', 'Reactor 3')

% plot control input
figure
grid on
xlabel('sim_step')
plot(u_log(1,:), 'linewidth', 1.5)
hold on
plot(u_log(2,:), 'linewidth', 1.5)
hold on
plot(u_log(3,:), 'linewidth', 1.5)
ylabel('Control input (J/s)')
legend('u1', 'u2', 'u3')

%% local functions

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
%global B_terminal_cost
J = 0; % initialize cost
if ORDER == 1
    u_overall = [u_host u_adj];
elseif ORDER == 2
    %u_overall = [u_adj(:,1) u_host u_adj(:,2)];
    u_overall = [u_host u_adj];
else
    u_overall = [u_adj u_host];
end

for i = 1:Np
    u_i = u_overall(min(i, Nc), :)';
    x_sample = sysfwd(subsys, x0, u_i); % size: size(x0, 2) x round(sample_interval/simulation_interval)
    % calculate cost in one sample interval
    delta_subx_array = x_sample - repmat(xs, 1, size(x_sample, 2)); % X - Xs
    sub_weighted_norm = trace(delta_subx_array' * Q * delta_subx_array);
    sub_discreted_normed_integeral = simulation_interval * sub_weighted_norm;
    J = J + sub_discreted_normed_integeral + u_i'*Ru*u_i* sample_interval;
    x0 = x_sample(:, end); % update state
end
end


% nonlinear constraint function
function [cineq, ceq] = nonlinear_horizon(subsys, u_host, u_adj, x0, N, ORDER)
% suppose ||x0-x_barrier||>=d_barrier
global B_terminal_cost Nc
if ~exist('N', 'var')
    N = 1;
end

% if ~exist('ORDER', 'var')
%     ORDER = 1;
% end

if ORDER == 1
    u_overall = [u_host u_adj];
elseif ORDER == 2
    %u_overall = [u_adj(:,1) u_host u_adj(:,2)];
    u_overall = [u_host u_adj];
else
    u_overall = [u_adj u_host];
end

cineq = zeros(N, 1); ceq = [];
for i = 1:N
    u_i = u_overall(min(i, Nc), :)';
    x0_sample = sysfwd(subsys, x0, u_i);
    x0_new = x0_sample(:, end); % new state
    cineq(i, :) = B_terminal_cost(x0_new);
    x0 = x0_new; % update state
end
end

function [cineq, ceq] = nonlinear_gamma_horizon(subsys, gamma, u_host, u_adj, u_host_pre, u_adj_pre, x0, N, ORDER)
% suppose ||x0-x_barrier||>=d_barrier
global B_terminal_cost Nc

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
    %u_overall = [u_adj(:,1) u_host u_adj(:,2)];
    u_overall = [u_host u_adj];
else
    u_overall = [u_adj u_host];
end

for i = 1:N
    u_i = u_overall(min(i, Nc), :)';
    x0_sample = sysfwd(subsys, x0, u_i);
    x0_new = x0_sample(:, end);
    cineq(i, :) = B_terminal_cost(x0_new);
    x0 = x0_new; % update state
end
end