clc
clear all
close all

load cmpc
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

load dmpc
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
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -250 400])
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
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -15 600])
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
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -50 600])
ylabel('u_3');
xlabel ('Time(h)')



com1=(dcost1-ccost1)/ccost1
com2=(dcost2-ccost2)/ccost2
com3=(dcost3-ccost3)/ccost3









