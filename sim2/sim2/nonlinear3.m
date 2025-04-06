function [cineq,ceq] = nonlinear3(u3,u3s,x7,x8,x9,x7s,x8s,x9s,x4,x5,x6,Delta,interval,pred,F2,F3,FD,FR,V3,alA,alB,alC,Cp)
x4s=6.41585647871797e-05;
x5s=0.000695133982828163;
x6s=475.482353426126;
   ceq = [];
%% adaptive controller
%    rho=0.5; 
%    alp42=0.8;  
%    alp52=0.6;
%    alp62=500;
%    alp72=0.62;
%    alp82=0.51;  
%    alp91=200;
%    alp92=500;
%    VX4=(x4-x4s)*(x4-x4s);
%    VX5=(x5-x5s)*(x5-x5s);
%    VX6=(x6-x6s)*(x6-x6s);
%    VX7=(x7-x7s)*(x7-x7s);
%    VX8=(x8-x8s)*(x8-x8s);
%    VX9=(x9-x9s)*(x9-x9s);
%    kx1=alp91/((1-rho)*alp91)/alp42*VX4;
%    kx2=alp91/((1-rho)*alp91)/alp52*VX5;
%    kx3=alp91/((1-rho)*alp91)/alp62*VX6;
%    kx4=alp91/((1-rho)*alp91)/alp72*VX7;
%    kx5=alp91/((1-rho)*alp91)/alp82*VX8;
%    kx6=alp91/((1-rho)*alp91)/alp92*VX9;
%    kx=max([kx1,kx2,kx3,kx4,kx5,kx6]);
%    
%    if(x9-x9s>0)
%        kxt=-kx;
%        kx=kxt;
%    end 
   
%% nonlinear
%    LfV=2*((x7-x7s)+(x8-x8s)+(x9-x9s))*((F2*x6-(FD+FR)*x9-F3*x9)/V3);
%    LgV=2*((x7-x7s)+(x8-x8s)+(x9-x9s))/(Cp*V3);
%    dVx=2*((x7-x7s)+(x8-x8s)+(x9-x9s))*((F2*x6-(FD+FR)*x9-F3*x9)/V3+kx/(Cp*V3));  
%    cineq = LfV+LgV*u3-dVx+0.1;
    alpha_3 = F3/V3;
    sigma = 2.806 * alpha_3;
    V11 = 0.0001*(x7-x7s)^2 + 0.0001*(x8-x8s)^2 + 0.5*(x9-x9s)^2;
    V12 = 0.0001*(x4-x4s)^2 + 0.0001*(x5-x5s)^2 + 0.5*(x6-x6s)^2;
    dV = 2 * [0.0001*(x7-x7s),0.0001*(x8-x8s),0.5*(x9-x9s)];
    Lf = dV * [(F2*x4-(FD+FR)*(alA*x7/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x7)/V3;(F2*x5-(FD+FR)*(alB*x8/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x8)/V3;(F2*x6-(FD+FR)*x9-F3*x9)/V3];
    Lg = dV * [0;0;1/(Cp*V3)];
%     cineq = Lf + Lg * u3 - alpha_3 * V;
    error = Lf + alpha_3 * V11 - sigma * V12
    cineq = Lf + Lg * u3 + alpha_3 * V11 - sigma * V12;
end

