function fun = SUB2(u2,u2s,x4,x5,x6,x4s,x5s,x6s,x1,x2,x3,Delta,interval,pred,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
fun = 0;

for j=1:pred
    
 for k = 1:Delta/interval
    
    eq2(1) = (Ff2*xA0+F1*x1-F2*x4)/V2-kA*exp(-EAR/x6)*x4;%x4
    eq2(2) = (F1*x2-F2*x5)/V2+kA*exp(-EAR/x6)*x4-kB*exp(-EBR/x6)*x5;%x5
    eq2(3) = (Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp+u2(j)/(Cp*V2); %x6
    

    x4 = x4 + interval*eq2(1);
    x5 = x5 + interval*eq2(2);
    x6 = x6 + interval*eq2(3);

    fun = fun + ((x4-x4s)*(x4-x4s)+(x5-x5s)*(x5-x5s)+(x6-x6s)*(x6-x6s)+0.0001*(u2(j)-u2s)*(u2(j)-u2s))*interval;

 end
 
end
