function fun = SUB(u,us,x,xs,Delta,interval,pred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
fun = 0;

for j=1:pred
    
 for k = 1:Delta/interval
    
    eq(1) = (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/x(3))*x(1);
    eq(2) = (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/x(3))*x(1)-kB*exp(-EBR/x(3))*x(2);
    eq(3) = (Ff1*T0+FR*x(9)-F1*x(3))/V1-(kA*exp(-EAR/x(3))*x(1)*dHA+kB*exp(-EBR/x(3))*x(2)*dHB)/Cp+u(j)/(Cp*V1);
    eq(4) = (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/x(6))*x(4);
    eq(5) = (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/x(6))*x(4)-kB*exp(-EBR/x(6))*x(5);
    eq(6) = (Ff2*T0+F1*x(3)-F2*x(6))/V2-(kA*exp(-EAR/x(6))*x(4)*dHA+kB*exp(-EBR/x(6))*x(5)*dHB)/Cp+u(j+pred)/(Cp*V2); 
    eq(7) = (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
    eq(8) = (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
    eq(9) = (F2*x(6)-(FD+FR)*x(9)-F3*x(9))/V3+u(j+2*pred)/(Cp*V3);
    
    x = x + interval*eq;

    %fun = fun + ((x(1)-xs(1))*(x(1)-xs(1))+(x(2)-xs(2))*(x(2)-xs(2))+10000*(x(3)-xs(3))*(x(3)-xs(3))+(x(4)-xs(4))*(x(4)-xs(4))+...
     %   (x(5)-xs(5))*(x(5)-xs(5))+10000*(x(6)-xs(6))*(x(6)-xs(6))+(x(7)-xs(7))*(x(7)-xs(7))+(x(8)-xs(8))*(x(8)-xs(8))+10000*(x(9)-xs(9))*(x(9)-xs(9))+...
      %  1000*(u(j)-us(1))*(u(j)-us(1))+1000*(u(j+pred)-us(2))*(u(j+pred)-us(2))+10000*(u(j+2*pred)-us(3))*(u(j+2*pred)-us(3)))*interval;
 
 fun = fun + (1.5*(x(1)-xs(1))*(x(1)-xs(1))+1.5*(x(2)-xs(2))*(x(2)-xs(2))+0.2*(x(3)-xs(3))*(x(3)-xs(3))+1.2*(x(4)-xs(4))*(x(4)-xs(4))+...
        1.2*(x(5)-xs(5))*(x(5)-xs(5))+0.18*(x(6)-xs(6))*(x(6)-xs(6))+1.15*(x(7)-xs(7))*(x(7)-xs(7))+1.15*(x(8)-xs(8))*(x(8)-xs(8))+0.15*(x(9)-xs(9))*(x(9)-xs(9))+...
        0.0001*(u(j)-us(1))*(u(j)-us(1))+0.0001*(u(j+pred)-us(2))*(u(j+pred)-us(2))+0.0001*(u(j+2*pred)-us(3))*(u(j+2*pred)-us(3)))*interval;
 end
 
end