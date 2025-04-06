function [cineq,ceq] = nonlinear2(u2,x1,x2,x3,x4,x5,x6,x1s,x2s,x3s,x4s,x5s,x6s,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0)

   ceq = [];
%% adaptive controller
   rho=0.5;
   alp12=0.83;
   alp22=0.65;
   alp32=500;
   alp42=0.8;  
   alp52=0.6;
   alp61=200;
   alp62=500;
   VX1=(x1-x1s)*(x1-x1s);
   VX2=(x2-x2s)*(x2-x2s);
   VX3=(x3-x3s)*(x3-x3s);
   VX4=(x4-x4s)*(x4-x4s);
   VX5=(x5-x5s)*(x5-x5s);
   VX6=(x6-x6s)*(x6-x6s);
   kx1=alp61/((1-rho)*alp61)/alp12*VX1;
   kx2=alp61/((1-rho)*alp61)/alp22*VX2;
   kx3=alp61/((1-rho)*alp61)/alp32*VX3;
   kx4=alp61/((1-rho)*alp61)/alp42*VX4;
   kx5=alp61/((1-rho)*alp61)/alp52*VX5;
   kx6=alp61/((1-rho)*alp61)/alp62*VX6;
   kx=max([kx1,kx2,kx3,kx4,kx5,kx6]);
   if(x6-x6s>0)
       kxt=-kx;
       kx=kxt;
   end 
   
%% nonlinear
   LfV=2*((x4-x4s)+(x5-x5s)+(x6-x6s))*((Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp);
   LgV=2*((x4-x4s)+(x5-x5s)+(x6-x6s))/(Cp*V2);
   dVx=2*((x4-x4s)+(x5-x5s)+(x6-x6s))*((Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp+kx/(Cp*V2));  
   cineq = LfV+LgV*u2-dVx+0.1;

end
