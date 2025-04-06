function [cineq,ceq] = nonlinear1(u1,u1s,x1,x2,x3,x1s,x2s,x3s,x7,x8,x9,x7s,x8s,x9s,Delta,interval,pred,u1pred,F1,FR,Ff1,V1,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)

   ceq = [];
%% adaptive controller
%    rho=0.5;
%    alp12=0.83;
%    alp22=0.65;
%    alp31=200;   
%    alp32=500;
%    alp72=0.62;
%    alp82=0.51;
%    alp92=500;
%    VX1=(x1-x1s)*(x1-x1s);
%    VX2=(x2-x2s)*(x2-x2s);
%    VX3=(x3-x3s)*(x3-x3s);
%    VX7=(x7-x7s)*(x7-x7s);
%    VX8=(x8-x8s)*(x8-x8s);
%    VX9=(x9-x9s)*(x9-x9s);
%    kx1=alp31/((1-rho)*alp31)/alp12*VX1;
%    kx2=alp31/((1-rho)*alp31)/alp22*VX2;
%    kx3=alp31/((1-rho)*alp31)/alp32*VX3;
%    kx4=alp31/((1-rho)*alp31)/alp72*VX7;
%    kx5=alp31/((1-rho)*alp31)/alp82*VX8;
%    kx6=alp31/((1-rho)*alp31)/alp92*VX9;
%    kx=max([kx1,kx2,kx3,kx4,kx5,kx6]);
% 
%    if(x3-x3s>0)
%        kxt=-kx;
%        kx=kxt;
%    end 
%    
%% nonlinear
%    LfV=2*((x1-x1s)+(x2-x2s)+(x3-x3s))*((Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp);
%    LgV=2*((x1-x1s)+(x2-x2s)+(x3-x3s))/(Cp*V1);
%    dVx=2*((x1-x1s)+(x2-x2s)+(x3-x3s))*((Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp)+2*((x1-x1s)+(x2-x2s)+(x3-x3s))*kx/(Cp*V1);  
%    cineq = LfV+LgV*u1-dVx+0.1;
    alpha_3 = F1/V1;
    sigma = 0.75 * alpha_3;
    V11 = 0.0001*(x1-x1s)^2 + 0.0001*(x2-x2s)^2 + 0.5*(x3-x3s)^2;
    V12 = 0.0001*(x7-x7s)^2 + 0.0001*(x8-x8s)^2 + 0.5*(x9-x9s)^2;
    dV = 2 * [0.0001*(x1-x1s),0.0001*(x2-x2s),0.5*(x3-x3s)];
    Lf = dV * [(Ff1*1+FR*(alA*x7/(alA*x7+alB*x8+alC*(1-x7-x8)))-F1*x1)/V1-kA*exp(-EAR/x3)*x1;(FR*(alB*x8/(alA*x7+alB*x8+alC*(1-x7-x8)))-F1*x2)/V1+kA*exp(-EAR/x3)*x1-kB*exp(-EBR/x3)*x2;(Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp];
    Lg = dV * [0;0;1/(Cp*V1)];
%     cineq = Lf + Lg * u1 + alpha_3 * V;
    error = Lf + alpha_3 * V11 - sigma * V12;
    cineq = Lf + Lg * u1 + alpha_3 * V11 - sigma * V12;
end

