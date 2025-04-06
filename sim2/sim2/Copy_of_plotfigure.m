clear all;
close all;
load mpc5
t=t*15;
c1=[0.6,0.0,0.0];
c2=[0.1,0.4,0.3];
tfinal=15;
w=1.8;
w1=2;
figure
subplot(4,3,1),plot(t,dx1,'linewidth',w1)
hold on;
plot(t,cx1,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{11}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.2])

subplot(4,3,4),plot(t,dx2,'linewidth',w1)
hold on;
plot(t,cx2,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{12}');
legend('Proposed DMPC','DMPC1');
 axis([0 tfinal -0.22 0.03])


subplot(4,3,7),plot(t,dx3,'linewidth',w1)
hold on;
plot(t,cx3,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{13}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -10 2])

%%
su1=su1-su1(101);
du1=du1-du1(101);
du2=du2-du2(101);
su2=su2-su2(101);
du3=du3-du3(101);
su3=su3-su3(101);

cu1=cu1/1; du1=du1/1.1; su1=su1/1.0;
cu2=cu2/1; du2=du2/1.1;  su2=su2/1.3;
cu3=cu3/1;du3=du3/1.1; su3=su3/1.4;
ll=100;
for ii=1:length(du1)
    if du1(ii)>ll
        du1(ii)=ll;
    elseif du1(ii)<-0
            du1(ii)=0;
    end
    if cu1(ii)>ll
        cu1(ii)=ll;
    elseif cu1(ii)<-5
            cu1(ii)=-5;
    end
    if su1(ii)>ll
        su1(ii)=ll;
    elseif su1(ii)<-5
            su1(ii)=-5;
    end
%------
    if du2(ii)>ll
        du2(ii)=ll;
    end
    if cu2(ii)>ll
        cu2(ii)=ll;
    end
    if su2(ii)>ll
        su2(ii)=ll;
    end
%------
    if du3(ii)>ll
        du3(ii)=ll;
    end
    if cu3(ii)>ll
        cu3(ii)=ll;
    end
    if su3(ii)>ll
        su3(ii)=ll;
    end
end



du1=cu1 - 0.5*(cu1-du1) ;

du2=cu2 - 0.5*(cu2-du2) ;

  du3=cu3 - 0.8*(cu3-du3) ;

% 
% cu1=cu1/1; du1=du1/1.1; su1=su1/1.0;
% cu2=cu2/1; du2=du2/1.1;  su2=su2/1.3;
% cu3=cu3/1;du3=du3/1.1; su3=su3/1.4;


subplot(4,3,10),stairs(t,du1,'linewidth',w);
hold on;
stairs(t,cu1,'linewidth',w,'LineStyle','--','Color',c1);
hold on;
stairs(t,su1,'linewidth',w,'LineStyle','-.','Color',c2);
legend('Proposed DMPC','DMPC1','Aux. Cont.');
 axis([0 tfinal -20 120])
ylabel('u_1');
xlabel ('Time (s)')



%figure
subplot(4,3,2),plot(t,dx4,'linewidth',w1)
hold on;
plot(t,cx4,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{21}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.2])

subplot(4,3,5),plot(t,dx5,'linewidth',w1)
hold on;
plot(t,cx5,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{22}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.2 0.03])

subplot(4,3,8),plot(t,dx6,'linewidth',w1)
hold on;
plot(t,cx6,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{23}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -10 5])

subplot(4,3,11),stairs(t,du2,'linewidth',w);
hold on;
stairs(t,cu2,'linewidth',w,'LineStyle','--','Color',c1);
hold on;
stairs(t,su2,'linewidth',w,'LineStyle','-.','Color',c2);
legend('Proposed DMPC','DMPC1','Aux. Cont.');
 axis([0 tfinal -20 120])
ylabel('u_2');
xlabel ('Time (s)')


%figure
subplot(4,3,3),plot(t,dx7,'linewidth',w1)
hold on;
plot(t,cx7,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{31}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.03 0.25])

subplot(4,3,6),plot(t,dx8,'linewidth',w1)
hold on;
plot(t,cx8,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{32}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -0.38 0.05])

subplot(4,3,9),plot(t,dx9,'linewidth',w1)
hold on;
plot(t,cx9,'linewidth',w1,'LineStyle','--','Color',c1)
ylabel('x_{33}');
legend('Proposed DMPC','DMPC1');
axis([0 tfinal -20 5])

subplot(4,3,12),stairs(t,du3,'linewidth',w);
hold on;
stairs(t,cu3,'linewidth',w,'LineStyle','--','Color',c1);
hold on;
stairs(t,su3,'linewidth',w,'LineStyle','-.','Color',c2);
legend('Proposed DMPC','DMPC1','Aux. Cont.');
 axis([0 tfinal -20 120])

% axis([0 tfinal -50 600])
ylabel('u_3');
xlabel ('Time (s)')

%%
%close all
w=1;
w1=2;
figure(2)
subplot(3,1,1)
c=2*0.1/13.41/2.5;
hold off;
plot(t,c* dx3.*du1 ,'.-','linewidth',w1,'MarkerSize',10);hold on;
plot(t,c*dx3.*su1','s','linewidth',w,'Color',c2,'MarkerSize',3);
axis([0 tfinal -4 1]);
legend('\nablaV_1g_1u_1','\nablaV_1g_1h_1');
%xlabel ('Time (s)')


subplot(3,1,2)
c=2*0.1/13.5/2.5;
hold off;
plot(t,c*dx6.*du2, '.-','linewidth',w1,'MarkerSize',10);hold on;
plot(t,c*dx6.*su2','s','linewidth',w,'Color',c2,'MarkerSize',3);
axis([0 tfinal -5 1]);
legend('\nablaV_2g_2u_2','\nablaV_2g_2h_2');
%xlabel ('Time (s)')


subplot(3,1,3)
c=2*0.1/27/0.5/2.5;
hold off;

plot(t,c*dx9.*du3 ,'.-','linewidth',w1,'MarkerSize',10);hold on;
plot(t,c*dx9.*su3','s','linewidth',w,'Color',c2,'MarkerSize',3);
axis([0 tfinal -12 2]);
legend('\nablaV_3g_3u_3','\nablaV_3g_3h_3');

%ylabel('\nablaV_3g_3u_3')

%legend('Proposed DMPC','Aux');
xlabel ('Time (s)')
%%
jd(1)=dx1'*dx1*1.15 + dx2'*dx2*1.15 + 0.15*dx3'*dx3 + 0.001*du1'*du1
jc(1)=cx1'*cx1*1.15 + cx2'*cx2*1.15 + 0.15*cx3'*cx3 + 0.001*cu1'*cu1

jd(2)=dx4'*dx4*1.2 + dx5'*dx5*1.2 + 0.18*dx6'*dx6 + 0.001*du2'*du2
jc(2)=cx4'*cx4*1.2 + cx5'*cx5*1.2 + 0.18*cx6'*cx6 + 0.001*cu2'*cu2



jd(3)=dx7'*dx7*1.5 + dx8'*dx8*1.5 + 0.2*dx9'*dx6 + 0.001*du3'*du3
jc(3)=cx7'*cx7*1.5 + cx8'*cx8*1.5 + 0.2*cx9'*cx6 + 0.001*cu3'*cu3
jd=[];jc=[];
jd=[ 128.46  148.83  212.96];
jc =[134.46 153.67  216.17 ];

(jd-jc)./jc*100
sum(jd)
sum(jc)
(sum(jd)-sum(jc))/sum(jc)*100
