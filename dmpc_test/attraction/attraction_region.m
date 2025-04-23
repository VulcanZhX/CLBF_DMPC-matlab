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

func1 = @(x, x_adj, u)[
    (Ff1*xA0+FR*(alA*x_adj(4)/(alA*x_adj(4)+alB*x_adj(5)+alC*(1-x_adj(4)-x_adj(5))))-F1*x(1))/V1-kA*exp(-EAR/(x(3)*1))*x(1);
    (FR*(alB*x_adj(5)/(alA*x_adj(4)+alB*x_adj(5)+alC*(1-x_adj(4)-x_adj(5))))-F1*x(2))/V1+kA*exp(-EAR/(x(3)*1))*x(1)-kB*exp(-EBR/(x(3)*1))*x(2);
    ((Ff1*T0+FR*(x_adj(6)*1)-F1*(x(3)*1))/V1-(kA*exp(-EAR/(x(3)*1))*x(1)*dHA+kB*exp(-EBR/(x(3)*1))*x(2)*dHB)/Cp+u(1)/(Cp*V1));
];

% solve system equilibrium
u_stable = [0 0 0]';
x0_initial = [0.33 0.65 640 0.265 0.61 521 0.32 0.670 445]'; x0 = x0_initial;
[xs, fval] = fsolve(@(x) func(x, u_stable), x0_initial, optimset('Display', 'off'));

% Attraction Region of S1
xs1 = xs(1:3);
A11_linear = [-(F1/V1+1.25*kA) 0 0;
            1.25*kA -(F1/V1+1.35*kB) 0;
            -1.25*kA*dHA/Cp -1.35*kB*dHB/Cp -F1/V1];

Q_lyap1 = eye(3);
P_lyap1 = lyap(A11_linear', Q_lyap1);
c1_sys1 = min(eig(P_lyap1));
c2_sys1 = max(eig(P_lyap1));

Vlyap = @(x)(x-xs1)'*P_lyap1*(x-xs1);
Vlyap_dt = @(x, x_adj, u)2*(x-xs1)'*P_lyap1*func1(x, x_adj, u);

%Vlyap(x0_initial(1:3))
%Vlyap_dt(x0_initial(1:3), x0_initial(4:end), u_stable)

opt_options = optimoptions('fmincon', 'Display', 'final-detailed', 'Algorithm', 'sqp');

nonlin_constraint = @(x)nonlin_con(@(x) Vlyap_dt(x, x0_initial(4:end), u_stable), x, xs1);
opt_problem_attraction1 = createOptimProblem('fmincon', 'objective', @(x) Vlyap(x), ...
    'x0', [0.12 0.23 400]', 'lb', [0 0 300]', 'ub', [1 1 800]', 'nonlcon', nonlin_constraint, 'options', opt_options);
[x_opt, V_opt] = fmincon(opt_problem_attraction1);

function [cineq, ceq] = nonlin_con(nonlinfunc, x, xs)
    % Nonlinear inequality constraints
    ceq = nonlinfunc(x);
    % Nonlinear equality constraints
    cineq = [0.02-norm(x(1:2)-xs(1:2)); 1-abs(xs(3)-x(3))];
end

Vlyap_dt([0.5 0.5 453.510]', x0_initial(4:end), u_stable)