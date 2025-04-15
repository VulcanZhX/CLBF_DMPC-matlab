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

% solve system equilibrium
u_stable = [0 0 0]';

xs1=0;
xs2=0;
xs3=499.479234575489;
xs4=0;
xs5=0;
xs6=475.482353426126;
xs7=0;
xs8=0;
xs9=314.757389283257;
x_initial  = [xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9]';
[x_eq, fval] = fsolve(@(x) func(x, u_stable), x_initial, optimset('Display', 'off'));
% show equilibrium state
fprintf('Equilibrium state:\n');
disp(norm(fval))