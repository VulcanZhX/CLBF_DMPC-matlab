function [cineq,ceq] = nonlinear1(u1,x1,x2,x3,x7,x8,x9,x1s,x2s,x3s,x7s,x8s,x9s,F1,FR,Ff1,V1,kA,kB,EAR,EBR,dHA,dHB,Cp,T0)

   ceq = [];
%% adaptive controller
   rho=0.5;
   alp12=0.83;
   alp22=0.65;
   alp31=200;   
   alp32=500;
   alp72=0.62;
   alp82=0.51;
   alp92=500;
   VX1=(x1-x1s)*(x1-x1s);
   VX2=(x2-x2s)*(x2-x2s);
   VX3=(x3-x3s)*(x3-x3s);
   VX7=(x7-x7s)*(x7-x7s);
   VX8=(x8-x8s)*(x8-x8s);
   VX9=(x9-x9s)*(x9-x9s);
   kx1=alp31/((1-rho)*alp31)/alp12*VX1;
   kx2=alp31/((1-rho)*alp31)/alp22*VX2;
   kx3=alp31/((1-rho)*alp31)/alp32*VX3;
   kx4=alp31/((1-rho)*alp31)/alp72*VX7;
   kx5=alp31/((1-rho)*alp31)/alp82*VX8;
   kx6=alp31/((1-rho)*alp31)/alp92*VX9;
   kx=max([kx1,kx2,kx3,kx4,kx5,kx6]);


   if(x3-x3s>0)
       kxt=-kx;
       kx=kxt;
   end 
   
%% nonlinear
   LfV=2*((x1-x1s)+(x2-x2s)+(x3-x3s))*((Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp);
   LgV=2*((x1-x1s)+(x2-x2s)+(x3-x3s))/(Cp*V1);
   dVx=2*((x1-x1s)+(x2-x2s)+(x3-x3s))*((Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp)+2*((x1-x1s)+(x2-x2s)+(x3-x3s))*kx/(Cp*V1);  
   cineq = LfV+LgV*u1-dVx+0.1;
   ceq = -0.023*x^2 + 0.23*abs(x_2)
end

