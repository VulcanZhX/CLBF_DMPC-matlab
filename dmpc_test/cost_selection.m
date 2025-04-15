clear
rng('default')

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
x0_initial = [0.23 0.35 590 0.42 0.47 450 0.31 0.68 341]'; x0 = x0_initial;

% time parameters
global sample_interval simulation_interval
sample_interval = 0.1; simulation_interval = 1e-4;
sim_steps = 20000; % number of simulation steps, overall time = sim_steps*sample_interval
max_p_iteration = 10; % number of iterations within one sample interval

% controller parameters
global Nc Np Q Ru
Nc = 1; % controller horizon
Np = 2; % prediction horizon
Q = diag([1, 1, 2, 1, 1, 2, 1, 1, 2]); % state cost
Ru = 1e-10*eye(3); % control cost

% initial input guess and bounds
U0 = 12.6e5 * ones(Nc, 3); % initial guess
U_lb = repmat([12.6e5-0.56e6 16.2e5-0.55e5 -12.6e5+0.55e5], Nc, 1); % lower bound
U_ub = repmat([12.6e5+0.56e6 16.2e5+0.55e5 12.6e5+0.55e5], Nc, 1); % upper bound

% log the state and control input
x_log = []; u_log = [];


x_nxt = sysfwd(func, x0_initial, [0 0 0]');
disp(x_nxt(:, end))

function x_nxt = sysfwd(sys, x0, u)
% sysfwd: forward simulation of the system (one sample interval)
global sample_interval simulation_interval
x_nxt = zeros(size(x0, 1), round(sample_interval/simulation_interval));
fwd_step = round(sample_interval/simulation_interval);
for i = 1:500*fwd_step
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
    J = J + sub_discreted_normed_integeral + u_i' * Ru * u_i * simulation_interval;
    x0 = x_sample(:, end); % update state
end
end