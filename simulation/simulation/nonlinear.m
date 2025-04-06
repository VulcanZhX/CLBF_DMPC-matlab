function [cineq,ceq] = nonlinear(u,x,xs,pred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)

   ceq = [];
%% adaptive controller
    alp1=5e-7;
    alp2=5e-7;
    alp3=5e-4;
    alp4=5e-7;
    alp5=5e-7;
    alp6=5e-4;
    alp7=5e-8;
    alp8=5e-8;
    alp9=5e-5;
    
    f1 = (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/x(3))*x(1);
    f2 = (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/x(3))*x(1)-kB*exp(-EBR/x(3))*x(2);
    f3 = (Ff1*T0+FR*x(9)-F1*x(3))/V1-(kA*exp(-EAR/x(3))*x(1)*dHA+kB*exp(-EBR/x(3))*x(2)*dHB)/Cp;
    f4 = (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/x(6))*x(4);
    f5 = (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/x(6))*x(4)-kB*exp(-EBR/x(6))*x(5);
    f6 = (Ff2*T0+F1*x(3)-F2*x(6))/V2-(kA*exp(-EAR/x(6))*x(4)*dHA+kB*exp(-EBR/x(6))*x(5)*dHB)/Cp; 
    f7 = (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    f8 = (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    f9 = (F2*x(6)-(FD+FR)*x(9)-F3*x(9))/V3;
    
    g3=1/(Cp*V1);
    g6=1/(Cp*V2);
    g9=1/(Cp*V3);

   VX1=(x(1)-xs(1))*(x(1)-xs(1));
   VX2=(x(2)-xs(2))*(x(2)-xs(2));
   VX3=(x(3)-xs(3))*(x(3)-xs(3));
   VX4=(x(4)-xs(4))*(x(4)-xs(4));
   VX5=(x(5)-xs(5))*(x(5)-xs(5));
   VX6=(x(6)-xs(6))*(x(6)-xs(6));
   VX7=(x(7)-xs(7))*(x(7)-xs(7));
   VX8=(x(8)-xs(8))*(x(8)-xs(8));
   VX9=(x(9)-xs(9))*(x(9)-xs(9));
   V=alp1*VX1+alp2*VX2+alp3*VX3+alp4*VX4+alp5*VX5+alp6*VX6+alp7*VX7+alp8*VX8+alp9*VX9;
   
   Lfv=2*alp1*(x(1)-xs(1))*f1+2*alp2*(x(2)-xs(2))*f2+2*alp3*(x(3)-xs(3))*f3+2*alp4*(x(4)-xs(4))*f4+2*alp5*(x(5)-xs(5))*f5+...
       2*alp6*(x(6)-xs(6))*f6+2*alp7*(x(7)-xs(7))*f7+2*alp8*(x(8)-xs(8))*f8+2*alp9*(x(9)-xs(9))*f9;
   Lgv3=2*alp3*(x(3)-xs(3))*g3;
   Lgv6=2*alp6*(x(6)-xs(6))*g6;
   Lgv9=2*alp9*(x(9)-xs(9))*g9;
   
   if(Lgv3~=0)
       hx3=-(Lgv3+sqrt((Lfv^2+Lgv3^4)))/Lgv3;
   else
       hx3=0;
   end
   
   if(hx3>500)
       hx3=100;
   end
   if(hx3<-500)
       hx3=-100;
   end
   
   
   
   if(Lgv6~=0)
       hx6=-(Lgv6+sqrt((Lfv^2+Lgv6^4)))/Lgv6;
   else
       hx6=0;
   end
   
   if(hx6>500)
       hx6=500;
   end
   if(hx6<-500)
       hx6=-500;
   end
   
   
   if(Lgv9~=0)
       hx9=-(Lgv9+sqrt((Lfv^2+Lgv9^4)))/Lgv9;
   else
       hx9=0;
   end
   
   if(hx9>500)
       hx9=500;
   end
   if(hx9<-500)
       hx9=-500;
   end
   
   
   
%% nonlinear
   dVx=Lgv3*hx3+Lgv6*hx6+Lgv9*hx9;
   cineq=Lgv3*u(1)+Lgv6*u(1+pred)+Lgv9*u(1+2*pred)-dVx;

end

