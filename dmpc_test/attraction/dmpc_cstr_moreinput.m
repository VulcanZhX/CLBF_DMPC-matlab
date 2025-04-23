clear
rng('default')

% system parameters
global F1 F2 F3 FD FR V1 V2 V3 alA alB alC kA kB EAR EBR dHA dHB Cp T0 xA0
F1=74.53; F2=75.03; F3=8.168;
FD=0.662; FR=66.2;
% Ff1=8.33; Ff2=0.5;
V1=13.41; V2=13.5; V3=0.5;
alA=3.5; alB=1.1;alC=0.5;
kA=0.2; kB=0.018;
EAR=100; EBR=150;
dHA=-40; dHB=-50;
Cp=2.5; T0=313; xA0=1;

% system formulation
% x = [xa1 xb1 T1 xa2 xb2 T2 xa3 xb3 T3]
% u = [Q1 Q2 Q3 Ff1 Ff2]

func = @(x, u)[
    (u(4)*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/(x(3)*1))*x(1);
    (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/(x(3)*1))*x(1)-kB*exp(-EBR/(x(3)*1))*x(2);
    ((u(4)*T0+FR*(x(9)*1)-F1*(x(3)*1))/V1-(kA*exp(-EAR/(x(3)*1))*x(1)*dHA+kB*exp(-EBR/(x(3)*1))*x(2)*dHB)/Cp+u(1)/(Cp*V1));
    (u(5)*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/(x(6)*1))*x(4);
    (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/(x(6)*1))*x(4)-kB*exp(-EBR/(x(6)*1))*x(5);
    (u(5)*T0+F1*(x(3)*1)-F2*(x(6)*1))/V2-(kA*exp(-EAR/(x(6)*1))*x(4)*dHA+kB*exp(-EBR/(x(6)*1))*x(5)*dHB)/Cp+u(2)/(Cp*V2);
    (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    (F2*(x(6)*1)-(FD+FR)*(x(9)*1)-F3*(x(9)*1))/V3+u(3)/(Cp*V3);
    ];

% solve system equilibrium
u_stable = [0 0 0 8.33 0.5]';
x0_init = [0.6 0.25 450 0.4 0.52 420 0.51 0.35 410]'; x0 = x0_init;
[xs, fval] = fsolve(@(x) func(x, u_stable), x0_init, optimset('Display', 'off'));


% time parameters
global sample_interval simulation_interval
sample_interval = 0.025; simulation_interval = 1e-4;
sim_steps = 500; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 10; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 3; % controller horizon
Np = 6; % prediction horizon
Q = diag([1, 1, 0.005, 1, 1, 0.005, 1.2, 1.2, 0.01]); % state cost
Ru = 1e-4*diag([0.1 0.1 0.1 0.5 0.5]); % control cost

% initialize control input
U0 = rand(Nc, 5); % initial guess for control input
U_lb = repmat([-500 -500 -600 2.5 0], Nc, 1); % lower bound
U_ub = repmat([800 800 1000 15 1.25], Nc, 1); % upper bound

% log the state and control input
x_log = []; u_log = [];

% barriers
x_barrier = [0.4 0.335]'; % barrier point (for reactor 2)
d_barrier = 0.085; % distance to barrier
global B_terminal_cost
B_terminal_cost = @(x) (d_barrier - norm(x([1, 2])-x_barrier));

optim_options = optimoptions('fmincon', 'Display', 'off');
simulation_timer = 0;
%% rolling optimization
for i_sim = 1:sim_steps
    tic;
    for i_p = 1:max_p_iteration
        % optimization problem for input u1, u2, u3
        U0_host = U0(:, 1:3);
        U0_adj = U0(:, 4:5);
        localcost_q = @(U) sub_cost_function(func, U, U0_adj, x0, xs, u_stable, 1);
        nonlin_con = @(U) nonlinear_horizon(func, U, U0_adj, x0, Np, 1);
        opt_problem1 = createOptimProblem('fmincon', 'objective', localcost_q, ...
            'x0', U0_host, 'lb', U_lb(:, 1:3), 'ub', U_ub(:, 1:3), 'nonlcon', nonlin_con, 'options', optim_options);
        [U_opt1, J_opt1] = fmincon(opt_problem1);

        % optimization problem for input Ff1, Ff2(u4, u5)
        U0_host = U0(:, 4:5);
        U0_adj = U0(:, 1:3);
        localcost_q = @(U) sub_cost_function(func, U, U0_adj, x0, xs, u_stable, 2);
        nonlin_con = @(U) nonlinear_horizon(func, U, U0_adj, x0, Np, 2);
        opt_problem2 = createOptimProblem('fmincon', 'objective', localcost_q, ...
            'x0', U0_host, 'lb', U_lb(:, 4:5), 'ub', U_ub(:, 4:5), 'nonlcon', nonlin_con, 'options', optim_options);
        [U_opt2, J_opt2] = fmincon(opt_problem2);

        % convex combination of two optimization problems
        gamma_init = [0.5 0.5];
        weightcost_obj = @(gamma) sub_cost_function(func, gamma(1)*U_opt1 + (1-gamma(1))*U0(:, 1:3), ...
                        gamma(2)*U_opt2 + (1-gamma(2))*U0(:, 4:5), x0, xs, u_stable, 1);
        nonlin_con_cmb = @(gamma) nonlinear_gamma_horizon(func, gamma, U_opt1, U_opt2, U0(:, 1:3), U0(:, 4:5), x0, Np, 1);
        opt_problem_cmb = createOptimProblem('fmincon', 'objective', weightcost_obj, ...
            'x0', gamma_init, 'lb', [0 0], 'ub', [1 1], 'nonlcon', nonlin_con_cmb, 'options', optim_options);
        [gamma_opt, J_opt_cmb] = fmincon(opt_problem_cmb);
        % update control input
        U_opt1 = gamma_opt(1)*U_opt1 + (1-gamma_opt(1))*U0(:, 1:3);
        U_opt2 = gamma_opt(2)*U_opt2 + (1-gamma_opt(2))*U0(:, 4:5);
        U0 = [U_opt1 U_opt2];
        % iteration terminal condition
        if norm(U_opt1-U0(:, 1:3)) < 1 && norm(U_opt2-U0(:, 4:5)) < 1e-2
            break;
        end
    end
    % update state
    [~, x_full_sol] = ode45(@(t, x) func(x, U0(1, :)'), [0 sample_interval], x0);
    x0_new = x_full_sol(end, :)';
    x_log = [x_log x0_new];
    u_log = [u_log U0(1, :)'];
    x0 = x0_new; % update state
    % display elapsed time and current error
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f, safety: %.4f, J_cmb: %.4f\n', ...
        i_sim, toc, norm(x0 - xs), B_terminal_cost(x0), J_opt_cmb);
    simulation_timer = simulation_timer + sample_interval;
    if norm(x0([1:2, 4:5, 7:8]) - xs([1:2, 4:5, 7:8])) < 0.02 && norm(x0([3, 6, 9])-xs([3, 6, 9])) < 2
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
end


%% plot results
x_log = [x0_init x_log]; % add initial state to log

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

% plot inputs
% plot u1-u3
figure
subplot(2, 1, 1)
plot(u_log(1, :), 'linewidth', 1.5)
hold on
plot(u_log(2, :), 'linewidth', 1.5)
hold on
plot(u_log(3, :), 'linewidth', 1.5)
ylabel('Control input')
legend('Q1', 'Q2', 'Q3')

% plot u4-u5
subplot(2, 1, 2)
plot(u_log(4, :), 'linewidth', 1.5)
hold on
plot(u_log(5, :), 'linewidth', 1.5)
xlabel('sim step')
ylabel('Control input')
legend('Ff1', 'Ff2')

%% Simulating autonomous evolution

x0 = x0_init;
[t, x_global_zero] = ode45(@(t, x)func(x, u_stable), [0 sim_steps*sample_interval], x0);
% interpolate the state
x_glboal_interp_zero = interp1(t, x_global_zero, 0:sample_interval:sim_steps*sample_interval);

% plot natural state evolution
figure
plot(x_glboal_interp_zero(:, 1), x_glboal_interp_zero(:, 2), 'linewidth', 1, 'Color', 'r')
hold on
plot(x_glboal_interp_zero(:, 4), x_glboal_interp_zero(:, 5), 'linewidth', 1, 'Color', 'g')
hold on
plot(x_glboal_interp_zero(:, 7), x_glboal_interp_zero(:, 8), 'linewidth', 1, 'Color', 'b')
hold on

scatter(xs(1), xs(2), 'r*')
hold on
scatter(xs(4), xs(5), 'g*')
hold on
scatter(xs(7), xs(8), 'b*')
hold on

% plot temperature
figure
grid on
xlabel('sim_step')
plot(x_glboal_interp_zero(:, 3), 'linewidth', 1.5)
hold on
plot(x_glboal_interp_zero(:, 6), 'linewidth', 1.5)
hold on
plot(x_glboal_interp_zero(:, 9), 'linewidth', 1.5)
ylabel('Temperature (K)')
legend('T1', 'T2', 'T3')

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

function J = sub_cost_function(subsys, u_host, u_adj, x0, xs, us, ORDER)
global sample_interval simulation_interval Nc Np Q Ru
% global B_terminal_cost
J = 0; % initialize cost

if ORDER == 1
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
    J = J + sub_discreted_normed_integeral + (u_i-us)'*Ru*(u_i-us)* sample_interval;
    x0 = x_sample(:, end); % update state
end
end

% nonlinear constraint function
function [cineq, ceq] = nonlinear_horizon(subsys, u_host, u_adj, x0, N, ORDER)
global B_terminal_cost Nc
if ~exist('N', 'var')
    N = 1;
end

if ~exist('ORDER', 'var')
    ORDER = 1;
end

if ORDER == 1
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
global B_terminal_cost Nc

if ~exist('N', 'var')
    N = 1;
end

if ~exist('ORDER', 'var')
    ORDER = 1;
end

cineq = zeros(N, 1); ceq = [];
u_host = gamma(1)*u_host + (1-gamma(1))*u_host_pre;
u_adj = gamma(2)*u_adj + (1-gamma(2))*u_adj_pre;
if ORDER == 1
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