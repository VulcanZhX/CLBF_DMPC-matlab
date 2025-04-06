clc
clear all
close all

format long
%% system parameters
Delta = 1e-2;
interval = 1e-4;
tfinal =1;
pred=5;


%% system parameters
F1=35;
F2=43;
F3=15;
FD=0.504;
FR=50.4;
Ff1=5.04;
Ff2=5.04;
V1=1*1000;
V2=0.5*1000;
V3=0.01*1000;
alA=3.5;
alB=1;
alC=0.5;
kA=2.77e3*5000;
kB=2.5e3*5000;
EA=50000;
EB=60000;
R=8.314;
EAR=EA/R;
EBR=EB/R;
MW=250e-3;
dHA=-60000/MW;
dHB=-70000/MW;
Cp=4.2e3;
T0=313;
xA0=1;

%% steady state
xs(1)=0.383;
xs(2)=0.581;
xs(3)=447.8;
xs(4)=0.391;
xs(5)=0.572;
xs(6)=444.6;
xs(7)=0.172;
xs(8)=0.748;
xs(9)=449.6;

xs1=2.74965498723184e-05;
xs2=0.000293813807189524;
xs3=499.479234575489;
xs4=6.41585647871797e-05;
xs5=0.000695133982828163;
xs6=475.482353426126;
xs7=7.41883444015044e-06;
xs8=0.000480500405435715;
xs9=314.757389283257;

xs = [xs1,xs2,xs3,xs4,xs5,xs6,xs7,xs8,xs9];



us(1)=0;
us(2)=0;
us(3)=0;

cost1=0;
cost2=0;
cost3=0;
%% initial state
x=[xs1+0.173;xs2-0.193;xs3-18;xs4+0.165;xs5-0.165;xs6-20;xs7+0.215;xs8-0.355;xs9-25];
% x=[0.225;0.110;176;0.187;0.053;105;0.063;0.110;135];


i = 1;
t(1) = 0;

%% state array
upred = zeros(3*pred,1);
uarray = zeros(tfinal/Delta+1,3*pred);


x1array = zeros(tfinal/Delta,1);
x2array = zeros(tfinal/Delta,1);
x3array = zeros(tfinal/Delta,1);
x4array = zeros(tfinal/Delta,1);
x5array = zeros(tfinal/Delta,1);
x6array = zeros(tfinal/Delta,1);
x7array = zeros(tfinal/Delta,1);
x8array = zeros(tfinal/Delta,1);
x9array = zeros(tfinal/Delta,1);

x1array(1,:) = x(1);
x2array(1,:) = x(2);
x3array(1,:) = x(3);
x4array(1,:) = x(4);
x5array(1,:) = x(5);
x6array(1,:) = x(6);
x7array(1,:) = x(7);
x8array(1,:) = x(8);
x9array(1,:) = x(9);




%% update
while i < (tfinal/Delta + 1)

    [u,upred,fval3] = controller(us,x,xs,Delta,interval,pred,upred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);
    uarray(i,:)=upred;
%       u(1)=0;
%       u(2)=0;
%       u(3)=0;
         
%% 真实系统
    for k = 1:Delta/interval
        eq(1) = (Ff1*xA0+FR*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(1))/V1-kA*exp(-EAR/x(3))*x(1);
        eq(2) = (FR*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F1*x(2))/V1+kA*exp(-EAR/x(3))*x(1)-kB*exp(-EBR/x(3))*x(2);
        eq(3) = (Ff1*T0+FR*x(9)-F1*x(3))/V1-(kA*exp(-EAR/x(3))*x(1)*dHA+kB*exp(-EBR/x(3))*x(2)*dHB)/Cp+u(1)/(Cp*V1);
        eq(4) = (Ff2*xA0+F1*x(1)-F2*x(4))/V2-kA*exp(-EAR/x(6))*x(4);
        eq(5) = (F1*x(2)-F2*x(5))/V2+kA*exp(-EAR/x(6))*x(4)-kB*exp(-EBR/x(6))*x(5);
        eq(6) = (Ff2*T0+F1*x(3)-F2*x(6))/V2-(kA*exp(-EAR/x(6))*x(4)*dHA+kB*exp(-EBR/x(6))*x(5)*dHB)/Cp+u(2)/(Cp*V2); 
        eq(7) = (F2*x(4)-(FD+FR)*(alA*x(7)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(7))/V3;
        eq(8) = (F2*x(5)-(FD+FR)*(alB*x(8)/(alA*x(7)+alB*x(8)+alC*(1-x(7)-x(8))))-F3*x(8))/V3;
        eq(9) = (F2*x(6)-(FD+FR)*x(9)-F3*x(9))/V3+u(3)/(Cp*V3);

        x(1) = x(1) + interval*eq(1);
        x(2) = x(2) + interval*eq(2);
        x(3) = x(3) + interval*eq(3);
        x(4) = x(4) + interval*eq(4);
        x(5) = x(5) + interval*eq(5);
        x(6) = x(6) + interval*eq(6);
        x(7) = x(7) + interval*eq(7);
        x(8) = x(8) + interval*eq(8);
        x(9) = x(9) + interval*eq(9);
    end
    
       
    x1array(i+1,:) =  x(1);
    x2array(i+1,:) =  x(2);
    x3array(i+1,:) =  x(3);
    x4array(i+1,:) =  x(4);
    x5array(i+1,:) =  x(5);
    x6array(i+1,:) =  x(6);
    x7array(i+1,:) =  x(7);
    x8array(i+1,:) =  x(8);
    x9array(i+1,:) =  x(9);

    t(i+1) = Delta* i;
    
    delx1=x(1)-xs(1);
    delx2=x(2)-xs(2);
    delx3=x(3)-xs(3);
    delx4=x(4)-xs(4);
    delx5=x(5)-xs(5);
    delx6=x(6)-xs(6);
    delx7=x(7)-xs(7);
    delx8=x(8)-xs(8);
    delx9=x(9)-xs(9);
    
    cost1=cost1+1.5*delx1*delx1+1.5*delx2*delx2+0.002*delx3*delx3+0.00001*u(1)*u(1)
    cost2=cost2+1.2*delx4*delx4+1.2*delx5*delx5+0.018*delx6*delx6+0.00001*u(2)*u(2)
    cost3=cost3+1.15*delx7*delx7+1.15*delx8*delx8+0.0015*delx9*delx9+0.00001*u(3)*u(3)


i = i +1
end

figure(1)
subplot(3,3,1),plot(t,x1array-xs1,'linewidth',1.5)
ylabel('x1');
subplot(3,3,2),plot(t,x2array-xs2,'linewidth',1.5)
ylabel('x2');
subplot(3,3,3),plot(t,x3array-xs3,'linewidth',1.5)
ylabel('x3');
subplot(3,3,4),plot(t,x4array-xs4,'linewidth',1.5)
ylabel('x4');
subplot(3,3,5),plot(t,x5array-xs5,'linewidth',1.5)
ylabel('x5');
subplot(3,3,6),plot(t,x6array-xs6,'linewidth',1.5)
ylabel('x6');
subplot(3,3,7),plot(t,x7array-xs7,'linewidth',1.5)
ylabel('x7');
subplot(3,3,8),plot(t,x8array-xs8,'linewidth',1.5)
ylabel('x8');
subplot(3,3,9),plot(t,x9array-xs9,'linewidth',1.5)
ylabel('x9');

xlabel ('Time(s)')
% title('Trajectory of x')

uarray(i,:) = upred;
figure(2)
subplot(3,1,1),stairs(t,uarray(:,1),'linewidth',1.5);
ylabel('u_1');
subplot(3,1,2),stairs(t,uarray(:,11),'linewidth',1.5);
ylabel('u_2');
subplot(3,1,3),stairs(t,uarray(:,21),'linewidth',1.5);
ylabel('u_3');
xlabel ('Time(s)')

save cmpc2;
