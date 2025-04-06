function [cineq,ceq] = nonlinear2(u2,u2s,x4,x5,x6,x4s,x5s,x6s,x1,x2,x3,Delta,interval,pred,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
    
x1s=2.74965498723184e-05;
x2s=0.000293813807189524;
x3s=499.479234575489;
   ceq = [];
%% adaptive controller
%    rho=0.5;
%    alp12=0.83;
%    alp22=0.65;
%    alp32=500;
%    alp42=0.8;  
%    alp52=0.6;
%    alp61=200;
%    alp62=500;
%    VX1=(x1-x1s)*(x1-x1s);
%    VX2=(x2-x2s)*(x2-x2s);
%    VX3=(x3-x3s)*(x3-x3s);
%    VX4=(x4-x4s)*(x4-x4s);
%    VX5=(x5-x5s)*(x5-x5s);
%    VX6=(x6-x6s)*(x6-x6s);
%    kx1=alp61/((1-rho)*alp61)/alp12*VX1;
%    kx2=alp61/((1-rho)*alp61)/alp22*VX2;
%    kx3=alp61/((1-rho)*alp61)/alp32*VX3;
%    kx4=alp61/((1-rho)*alp61)/alp42*VX4;
%    kx5=alp61/((1-rho)*alp61)/alp52*VX5;
%    kx6=alp61/((1-rho)*alp61)/alp62*VX6;
%    kx=max([kx1,kx2,kx3,kx4,kx5,kx6]);
%    if(x6-x6s>0)
%        kxt=-kx;
%        kx=kxt;
%    end 
   
%% nonlinear
%    LfV=2*((x4-x4s)+(x5-x5s)+(x6-x6s))*((Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp);
%    LgV=2*((x4-x4s)+(x5-x5s)+(x6-x6s))/(Cp*V2);
%    dVx=2*((x4-x4s)+(x5-x5s)+(x6-x6s))*((Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp+kx/(Cp*V2));  
%    cineq = LfV+LgV*u2-dVx+0.1;
    alpha_3 = F2/V2;
    sigma = 0.81 * alpha_3;
    V11 = 0.0001*(x4-x4s)^2 + 0.0001*(x5-x5s)^2 + 0.5*(x6-x6s)^2;
    V12 = 0.0001*(x1-x1s)^2 + 0.0001*(x2-x2s)^2 + 0.5*(x3-x3s)^2;
    dV = 2 * [0.0001*(x4-x4s),0.0001*(x5-x5s),0.5*(x6-x6s)];
    Lf = dV * [(Ff2*xA0+F1*x1-F2*x4)/V2-kA*exp(-EAR/x6)*x4;(F1*x2-F2*x5)/V2+kA*exp(-EAR/x6)*x4-kB*exp(-EBR/x6)*x5;(Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp];
    Lg = dV * [0;0;1/(Cp*V2)];
    cineq = Lf + Lg * u2 + alpha_3 * V11 - sigma * V12;
end
