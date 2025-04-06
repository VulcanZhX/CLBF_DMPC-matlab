clear all;
close all;
load mpc5
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

