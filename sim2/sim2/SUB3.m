function fun = SUB3(u3,u3s,x7,x8,x9,x7s,x8s,x9s,x4,x5,x6,Delta,interval,pred,F2,F3,FD,FR,V3,alA,alB,alC,Cp)
fun = 0;

for j=1:pred
    
 for k = 1:Delta/interval
    
    eq3(1) = (F2*x4-(FD+FR)*(alA*x7/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x7)/V3;
    eq3(2) = (F2*x5-(FD+FR)*(alB*x8/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x8)/V3;
    eq3(3) = (F2*x6-(FD+FR)*x9-F3*x9)/V3+u3(j)/(Cp*V3);
    
    x7 = x7 + interval*eq3(1);
    x8 = x8 + interval*eq3(2);
    x9 = x9 + interval*eq3(3);

    fun = fun + ((x7-x7s)*(x7-x7s)+(x8-x8s)*(x8-x8s)+(x9-x9s)*(x9-x9s)*1e2)*interval;

 end
 
end