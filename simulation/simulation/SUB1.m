function fun = SUB1(u1,u1s,x1,x2,x3,x1s,x2s,x3s,x7,x8,x9,Delta,interval,pred,F1,FR,Ff1,V1,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
fun = 0;

for j=1:pred
    
 for k = 1:Delta/interval
     % dx1=0.15;
     % dx2=500;         %有下标
    
    eq1(1) = (Ff1*xA0+FR*(alA*(x7)/(alA*(x7)+alB*(x8)+alC*(1-(x7)-(x8))))-F1*x1)/V1-kA*exp(-EAR/(x3))*x1;%x1
    eq1(2) = (FR*(alB*(x8)/(alA*(x7)+alB*(x8)+alC*(1-(x7)-(x8))))-F1*x2)/V1+kA*exp(-EAR/(x3))*(x1)-kB*exp(-EBR/(x3))*x2;%x2
    eq1(3) = (Ff1*T0+FR*(x9)-F1*x3)/V1-(kA*exp(-EAR/x3)*(x1)*dHA+kB*exp(-EBR/x3)*(x2)*dHB)/Cp+u1(j)/(Cp*V1);%x3
    
    x1 = x1 + interval*eq1(1);
    x2 = x2 + interval*eq1(2);
    x3 = x3 + interval*eq1(3);

    fun = fun + ((x1-x1s)*(x1-x1s)+(x2-x2s)*(x2-x2s)+(x3-x3s)*(x3-x3s)+0.0001*(u1(j)-u1s)*(u1(j)-u1s))*interval;

 end
 
end
