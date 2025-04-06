clc
clear all
close all

load cmpc5
xs1=2.74965498723184e-05;
xs2=0.000293813807189524;
xs3=490.485791311134;
xs4=6.41585647871797e-05;
xs5=0.000695133982828163;
xs6=463.932407617157;
xs7=7.41883444015044e-06;
xs8=0.000480500405435715;
xs9=307.008209455795;
cx1=x1array-xs1;
cx2=x2array-xs2;
cx3=x3array-xs3;
cx4=x4array-xs4;
cx5=x5array-xs5;
cx6=x6array-xs6;
cx7=x7array-xs7;
cx8=x8array-xs8;
cx9=x9array-xs9;
cu1=uarray(:,1);
cu2=uarray(:,11);
cu3=uarray(:,21);
ccost1=0;
ccost2=0;
ccost3=0;

load dmpc4
xs1=2.74965498723184e-05;
xs2=0.000293813807189524;
xs3=490.485791311134;
xs4=6.41585647871797e-05;
xs5=0.000695133982828163;
xs6=463.932407617157;
xs7=7.41883444015044e-06;
xs8=0.000480500405435715;
xs9=307.008209455795;
dx1=x1array-xs1;
dx2=x2array-xs2;
dx3=x3array-xs3;
dx4=x4array-xs4;
dx5=x5array-xs5;
dx6=x6array-xs6;
dx7=x7array-xs7;
dx8=x8array-xs8;
dx9=x9array-xs9;
du1=u1array(:,1);
du2=u2array(:,1);
du3=u3array(:,1);
us1 = du1(101);
us2 = du2(101);
us3 = du3(101);
us = [us1,us2,us3];
%     x1 = x1array(i);
%     x2 = x2array(i);
%     x3 = x3array(i);
%     x4 = x4array(i);
%     x5 = x5array(i);
%     x6 = x6array(i);
%     x7 = x7array(i);
%     x8 = x8array(i);
%     x9 = x9array(i);
%                                             
%     x = [x1,x2,x3,x4,x5,x6,x7,x8,x9]; 
%     xar = 3.5*x(7)/(3.5*x(7)+1*x(8)+0.5*(1-x(7)-x(8))); 
%     xbr = 1*x(8)/(3.5*x(7)+1*x(8)+0.5*(1-x(7)-x(8)));
%     partialV1=[2*1*(x(1)-xs(1)) 2*1*(x(2)-xs(2)) 2*0.0001*(x(3)-xs(3))];
%     Lf1 = partialV1*[(Ff1*xA0+FR*(alA*xar)-F1*x1)/V1-kA*exp(-EAR/x3)*x1;(FR*(xbr)-F1*x2)/V1+kA*exp(-EAR/x3)*x1-kB*exp(-EBR/x3)*x2;(Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp;];
%     Lg1 = partialV1*[0;0;1/(Cp*V1);];
%     partialV2=[2*1*(x(4)-xs(4)) 2*1*(x(5)-xs(5)) 2*0.0001*(x(6)-xs(6))];
%     Lf2 = partialV2*[(Ff2*xA0+F1*x1-F2*x4)/V2-kA*exp(-EAR/x6)*x4;(F1*x2-F2*x5)/V2+kA*exp(-EAR/x6)*x4-kB*exp(-EBR/x6)*x5;(Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp];
%     Lg2 = partialV2*[0;0;1/(Cp*V2);];
%     partialV3=[2*1*(x(7)-xs(7)) 2*1*(x(8)-xs(8)) 2*0.0001*(x(9)-xs(9))];
%     Lf3 = partialV3*[(F2*x4-(FD+FR)*(xar)-F3*x7)/V3;(F2*x5-(FD+FR)*(xbr)-F3*x8)/V3;(F2*x6-(FD+FR)*x9-F3*x9)/V3];
%     Lg3 = partialV3*[0;0;1/(Cp*V3);];    
%     su1(i) = ustd(Lf1,Lg1) + us(1);
%     su2(i) = ustd(Lf2,Lg2) + us(2);
%     su3(i) = ustd(Lf3,Lg3) + us(3);
dcost1=0;
dcost2=0;
dcost3=0;

i=1;

while i < 101
    ccost1=ccost1+1.15*cx1(i)*cx1(i)+1.15*cx2(i)*cx2(i)+0.0015*cx3(i)*cx3(i)+0.0000002*cu1(i)*cu1(i);
    ccost2=ccost2+1.2*cx4(i)*cx4(i)+1.2*cx5(i)*cx5(i)+0.018*cx6(i)*cx6(i)+0.0000007*cu2(i)*cu2(i);
    ccost3=ccost3+1.5*cx7(i)*cx7(i)+1.5*cx8(i)*cx8(i)+0.002*cx9(i)*cx9(i)+0.000003*cu3(i)*cu3(i);
    
    dcost1=dcost1+1.15*dx1(i)*dx1(i)+1.15*dx2(i)*dx2(i)+0.0015*dx3(i)*dx3(i)+0.0000002*du1(i)*du1(i);
    dcost2=dcost2+1.2*dx4(i)*dx4(i)+1.2*dx5(i)*dx5(i)+0.018*dx6(i)*dx6(i)+0.0000007*du2(i)*du2(i);
    dcost3=dcost3+1.5*dx7(i)*dx7(i)+1.5*dx8(i)*dx8(i)+0.002*dx9(i)*dx9(i)+0.000003*du3(i)*du3(i);
    
    i=i+1;
end



figure
subplot(4,1,1),plot(t,dx1,'linewidth',2)
hold on;
plot(t,cx1,'linewidth',2)
ylabel('x_{11}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.2])
subplot(4,1,2),plot(t,dx2,'linewidth',2)
hold on;
plot(t,cx2,'linewidth',2)
ylabel('x_{12}');
legend('Proposed DMPC','DMPC1');
 axis([0 tfinal -0.22 0.03])
subplot(4,1,3),plot(t,dx3,'linewidth',2)
hold on;
plot(t,cx3,'linewidth',2)
ylabel('x_{13}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -10 2])
subplot(4,1,4),stairs(t,du1,'linewidth',2);
hold on;
stairs(t,cu1,'linewidth',2);
hold on;
stairs(t,su1,'linewidth',2);
legend('Proposed DMPC','DMPC1','Aux');
% axis([0 tfinal -250 400])
ylabel('u_1');
xlabel ('Time(h)')

figure
subplot(4,1,1),plot(t,dx4,'linewidth',2)
hold on;
plot(t,cx4,'linewidth',2)
ylabel('x_{21}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.2])
subplot(4,1,2),plot(t,dx5,'linewidth',2)
hold on;
plot(t,cx5,'linewidth',2)
ylabel('x_{22}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.2 0.03])
subplot(4,1,3),plot(t,dx6,'linewidth',2)
hold on;
plot(t,cx6,'linewidth',2)
ylabel('x_{23}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -10 5])
subplot(4,1,4),stairs(t,du2,'linewidth',2);
hold on;
stairs(t,cu2,'linewidth',2);
hold on;
stairs(t,su2,'linewidth',2);
legend('Proposed DMPC','DMPC1','Aux');
% axis([0 tfinal -15 600])
ylabel('u_2');
xlabel ('Time(h)')

figure
subplot(4,1,1),plot(t,dx7,'linewidth',2)
hold on;
plot(t,cx7,'linewidth',2)
ylabel('x_{31}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.25])
subplot(4,1,2),plot(t,dx8,'linewidth',2)
hold on;
plot(t,cx8,'linewidth',2)
ylabel('x_{32}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.38 0.05])
subplot(4,1,3),plot(t,dx9,'linewidth',2)
hold on;
plot(t,cx9,'linewidth',2)
ylabel('x_{33}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -20 5])
subplot(4,1,4),stairs(t,du3,'linewidth',2);
hold on;
stairs(t,cu3,'linewidth',2);
hold on;
stairs(t,su3,'linewidth',2);
legend('Proposed DMPC','DMPC1','Aux');
% axis([0 tfinal -50 600])
ylabel('u_3');
xlabel ('Time(h)')



com1=(dcost1-ccost1)/ccost1
com2=(dcost2-ccost2)/ccost2
com3=(dcost3-ccost3)/ccost3
 save mpc4







