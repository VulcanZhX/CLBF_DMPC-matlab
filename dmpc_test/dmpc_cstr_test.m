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

% solve system equilibrium
u_stable = [0 0 0]';
x_initial  = [0.2 0.5 467 0.4 0.5 434 0.25 0.6 459]';
[x_eq, fval] = fsolve(@(x) func(x, u_stable), x_initial, optimset('Display', 'off'));
% show equilibrium state
fprintf('Equilibrium state:\n');
x_eq([3, 6, 9])
fval