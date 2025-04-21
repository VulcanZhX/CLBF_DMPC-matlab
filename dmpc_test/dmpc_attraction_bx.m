clear
rng('default')

%% Initialization

% system definition:
global func
func = @(x,u) [
    x(1) - x(1)^3 + x(2) + u(1);
    x(1) - 3*x(2) + u(2);
    ];
% func = @(x, u)[
%     - 3*x(1) + x(2) + u;
%     x(1) - 3*x(2)
%     ];

% initial state and target state
x0 = [-1.5 0.3]'; xs = [2/sqrt(3) 2/(sqrt(3)*3)]'; x_init = x0;

% Lyapunov function parameters
P = [3/16 1/16; 1/16 3/16];
A = [-3 1; 1 -3];
Lyap_dVt = @(x, u) (x-xs)'*P*func(x, u);

% time parameters
global sample_interval simulation_interval
sample_interval = 0.1; % sampling time
simulation_interval = 1e-3; % simulation time
sim_steps = 1000; % number of simulation steps, overall time = sim_steps*sample_interval
global Nc Np Q Ru
Q = eye(2); % state cost
Ru = 0.1; % input cost
Nc = 2; % controller horizon
Np = 4; % prediction horizon

U_initial = [0 0; 0 0];
U_lb = -3*ones(Nc, 2); % lower bound
U_ub = 3*ones(Nc, 2); % upper bound

% log the state and input
x_log = [];
u_log = [];

% control lyapunov-barrier function
d_barrier = 0.5; x_barrier = [-0.3 0]';

eps = 0.2; % abs(barrier_lowbound)
bar = @(x)max(d_barrier - sqrt((x-x_barrier)'*diag([1, 1])*(x-x_barrier)), -eps);
lyap_xs = @(x)0.25*(x-xs)'*P*(x-xs);

c1 = 0.25*min(eig(P));
c2 = 0.25*max(eig(P));
c3 = norm(x_barrier - xs) + d_barrier + eps;
c4 = norm(x_barrier - xs) - d_barrier;

% solve mu and kappa
mu = 3*(c2*c3-c1*c4)/eps;
kappa = -c1*c4;

% clbf definition
W_clbf = @(x)lyap_xs(x) + mu*bar(x) + kappa;

opt_option = optimoptions('fmincon', 'Display', 'off', 'Algorithm','interior-point','EnableFeasibilityMode',true);
simulation_timer = 0;
%% Rolling Optimization
for i_sim = 1:sim_steps
    tic;
    % define optimization problem

    % cost function
    J_cost = @(U) cost_function(U, x0, xs);
    % nonlinear constraints
    alpha_3x = 0.5; % Lyapunov function constraint
    nonlin= @(U) nonlinear_horizon(func, x0, U, @(x)W_clbf(x), Np);
    opt_problem = createOptimProblem('fmincon', 'objective', J_cost, ...
        'lb', U_lb, 'ub', U_ub, 'x0', U_initial, 'nonlcon', nonlin, 'options', opt_option);

    [U_opt, J_opt] = fmincon(opt_problem);

    % next state under U_opt(1)
    [t,x_sol]=ode45(@(t,x) func(x,U_opt(1, :)'), [0 sample_interval], x0);
    x0_new = x_sol(end,:)'; % new state
    x_log = [x_log x0]; % log the state
    u_log = [u_log U_opt(1,:)']; % log the control input
    % display elapsed time of the current loop
    fprintf('Simulation step %d, elapsed time: %.4f seconds, current error: %.4f, safety: %.4f\n',...
        i_sim, toc, norm(x0 - xs), W_clbf(x0));
    % terminal condition
    if norm(x0 - xs) < 0.01
        fprintf('Terminal condition met at step %d\n', i_sim);
        break;
    end
    simulation_timer = simulation_timer + sample_interval;
    x0 = x0_new; % update state
    U_initial = U_opt;   % update input
end

W_log = []; % log the Lyapunov function
for i = 1:size(x_log, 2)
    W_log = [W_log W_clbf(x_log(:, i))];
end


% plot input
% figure;
% plot(0:sample_interval:simulation_timer, u_log(1,:), 'g--', 'LineWidth', 1);
% grid on
% xlabel('Time (s)')
% ylabel('Control Input')
% title('Control Input over Time')

%% Auxiliary Stablizing Control
% presented by Sontag 1989-1991, auxiliary universal controller
func_f = @(x)[
    x(1) - x(1)^3 + x(2);
    x(1) - 3*x(2)
    ];

func_g = @(x)[1 0; 0 1];
LW = @(x)grad_Wx(W_clbf, x, 0.25);
LfW = @(x) LW(x)*func_f(x);
LgW = @(x) LW(x)*func_g(x);
gamma_gain = 6;
sontag_h = @(p,q,r) -(p+sqrt(p^2+r*norm(q)^4))*q'/(norm(q)^2);
sontag_combined = @(x) sontag_h(LfW(x), LgW(x), gamma_gain);
% sontag_combined(x0)

% Auxiliary Control Simulation
x_aux_log = []; u_aux_log = [];
x0_aux = x_init;
[t_aux, x_aux] = ode45(@(t, x)func(x, sontag_combined(x)), [0 sim_steps*sample_interval], x0_aux);
x_aux_fixed_interval = interp1(t_aux, x_aux, (0:simulation_interval:sim_steps*sample_interval-simulation_interval));


%% Auxiliary CLBF Control
% originated from Stabilization with guaranteed safety using Control
% Lyapunov–Barrier Function 2016. by Romdlony, etc.


%% Plot the results
figure;
% plot the trajectory
plot(x_log(1,:), x_log(2,:), 'r-', 'LineWidth', 2);
hold on; grid on
scatter(xs(1), xs(2), 'b*')
hold on
% plot the trajectory of the auxiliary system controller
plot(x_aux_fixed_interval(:,1), x_aux_fixed_interval(:,2), 'r-.', 'LineWidth', 1.5);
hold on
% plot lyapunov function derivative boundary
fimplicit(@(x1, x2) Lyap_dVt([x1; x2], [0 0]'), 'k--', 'LineWidth', 1.5);
fimplicit(@(x1, x2) bar([x1; x2]), 'g--', 'LineWidth', 1.5);
axis equal
legend('Trajectory', 'Target State', 'Stable Region Boundary')
xlabel('x_1')
ylabel('x_2')
title('Trajectory of the System')

% plot the Lyapunov function
figure;
plot(W_log, 'r-', 'LineWidth', 2);

%% Cost Function
function J = cost_function(U, x0, xs)
global func
global Nc Np Q Ru
global sample_interval simulation_interval

J = 0.0;
for i = 1:Np
    u_i = U(min(i, Nc), :)'; % get the control input
    [t,x_i]=ode45(@(t,x) func(x,u_i), [0 sample_interval], x0);
    x_i_fixed_interval = interp1(t, x_i, (0:simulation_interval:sample_interval-simulation_interval));
    % calculate cost in one sample interval
    delta_x_array = x_i_fixed_interval-repmat(xs', size(x_i_fixed_interval, 1), 1); %X - Xs
    discreted_normed_integeral = trace(delta_x_array*Q*delta_x_array')*simulation_interval;
    J = J + discreted_normed_integeral + u_i'*Ru*u_i*sample_interval;
    x0 = x_i_fixed_interval(end,:)';
end
% add terminal cost
end
%% Lyapunov_fun_cons
% dV/dt <= -a*V
function ineq_val = lyap_fun_cons(lyap_f, lyap_dvt, x, alpha)
    ineq_val = lyap_dvt(x) + alpha*lyap_f(x);
end

%% Derivative of 

%% Nonlinear Constriants
function [cineq, ceq] = nonlinear_horizon(subsys, x0, U, fun_con, N)
global sample_interval Nc
if ~exist('N', 'var')
    N = 1;
end
cineq = zeros(N, 1); ceq = [];

for i = 1:N
    u_i = U(min(i, Nc), :);
    cineq(i, :) = fun_con(x0);
    [~, x0_new] = ode45(@(t, x) subsys(x, u_i), [0 sample_interval], x0);
    x0_new = x0_new(end, :)';
    x0 = x0_new; % update state
end
end


function LW = grad_Wx(W_clbf, x, deltaxy)
    LW_x = W_clbf([x(1)+deltaxy x(2)]')-W_clbf(x);
    LW_y = W_clbf([x(1) x(2)+deltaxy]')-W_clbf(x);
    LW = [LW_x LW_y]/deltaxy;
end