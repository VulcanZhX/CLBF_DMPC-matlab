clc
clear all
close all

format long
%% system parameters
Delta = 1e-2;
interval = 1e-4;
tfinal =1;
pred=10;


%% system parameters
F1=35.5;
F2=43.5;
F3=15.5;
FD=0.504;
FR=50.4;
Ff1=5;
Ff2=5;
V1=1*1000;
V2=0.5*1000;
V3=0.012*1000;
alA=3.5;
alB=1;
alC=0.5;
kA=2.77e3*3600;
kB=2.5e3*3600;
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
x1s=0.383;
x2s=0.581;
x3s=447.8;
x4s=0.391;
x5s=0.572;
x6s=444.6;
x7s=0.172;
x8s=0.748;
x9s=449.6;

xs1=2.74965498723184e-05;
xs2=0.000293813807189524;
xs3=499.479234575489;
xs4=6.41585647871797e-05;
xs5=0.000695133982828163;
xs6=475.482353426126;
xs7=7.41883444015044e-06;
xs8=0.000480500405435715;
xs9=314.757389283257;


u1s=0;
u2s=0;
u3s=0;

cost1=0;
cost2=0;
cost3=0;
cost1_log = [];
cost2_log = [];
cost3_log = [];


%% initial state
% x1=0.529;
% x2=0.350;
% x3=325;
% x4=0.79;
% x5=0.36;
% x6=295;
% x7=0.10;
% x8=0.18;
% x9=276;

x1=xs1+0.173;
x2=xs2-0.193;
x3=xs3-18;
x4=xs4+0.165;
x5=xs5-0.165;
x6=xs6-20;
x7=xs7+0.215;
x8=xs8-0.355;
x9=xs9-25;
% x=[0.729;0.350;325;0.59;0.36;295;0.35;0.55;276];


i = 1;
t(1) = 0;

%% state array
u1pred = zeros(pred,1);
u2pred = zeros(pred,1);
u3pred = zeros(pred,1);
u1array = zeros(tfinal/Delta+1,pred);
u2array = zeros(tfinal/Delta+1,pred);
u3array = zeros(tfinal/Delta+1,pred);

x1array = zeros(tfinal/Delta,1);
x2array = zeros(tfinal/Delta,1);
x3array = zeros(tfinal/Delta,1);
x4array = zeros(tfinal/Delta,1);
x5array = zeros(tfinal/Delta,1);
x6array = zeros(tfinal/Delta,1);
x7array = zeros(tfinal/Delta,1);
x8array = zeros(tfinal/Delta,1);
x9array = zeros(tfinal/Delta,1);

x1array(1,:) = x1;
x2array(1,:) = x2;
x3array(1,:) = x3;
x4array(1,:) = x4;
x5array(1,:) = x5;
x6array(1,:) = x6;
x7array(1,:) = x7;
x8array(1,:) = x8;
x9array(1,:) = x9;




%% update
while i < (tfinal/Delta + 1)
    if(i>1)
        dx1=x1array(i-1,:); %第一个量 移位序列
        dx2=x2array(i-1,:);
        dx3=x3array(i-1,:);
        dx4=x4array(i-1,:);
        dx5=x5array(i-1,:);
        dx6=x6array(i-1,:);
        dx7=x7array(i-1,:);
        dx8=x8array(i-1,:);
        dx9=x9array(i-1,:);
    else
        dx1=x1; %sequence
        dx2=x2;
        dx3=x3;
        dx4=x4;
        dx5=x5;
        dx6=x6;
        dx7=x7;
        dx8=x8;
        dx9=x9;        
    end
    for iter=1:1

    [u1,u1pred,fval3] = CON1(u1s,x1,x2,x3,x1s,x2s,x3s,dx7,dx8,dx9,x7s,x8s,x9s,Delta,interval,pred,u1pred,F1,FR,Ff1,V1,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);
    u1array(i,:)=u1pred;
    
    [u2,u2pred,fval6] = CON2(u2s,x4,x5,x6,x4s,x5s,x6s,dx1,dx2,dx3,x1s,x2s,x3s,Delta,interval,pred,u2pred,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);
    u2array(i,:)=u2pred;
    
    [u3,u3pred,fval9] = CON3(u3s,x7,x8,x9,x7s,x8s,x9s,dx4,dx5,dx6,x4s,x5s,x6s,Delta,interval,pred,u3pred,F2,F3,FD,FR,V3,alA,alB,alC,Cp);
    u3array(i,:)=u3pred;
   
    % calcualte 
    % dx1-dx9 sequence

    end

    
    dx=10;
    

         
%% 真实系统
    for k = 1:Delta/interval
        eq(1) = (Ff1*xA0+FR*(alA*x7/(alA*x7+alB*x8+alC*(1-x7-x8)))-F1*x1)/V1-kA*exp(-EAR/x3)*x1;
        eq(2) = (FR*(alB*x8/(alA*x7+alB*x8+alC*(1-x7-x8)))-F1*x2)/V1+kA*exp(-EAR/x3)*x1-kB*exp(-EBR/x3)*x2;
        eq(3) = (Ff1*T0+FR*x9-F1*x3)/V1-(kA*exp(-EAR/x3)*x1*dHA+kB*exp(-EBR/x3)*x2*dHB)/Cp+u1/(Cp*V1);
        eq(4) = (Ff2*xA0+F1*x1-F2*x4)/V2-kA*exp(-EAR/x6)*x4;
        eq(5) = (F1*x2-F2*x5)/V2+kA*exp(-EAR/x6)*x4-kB*exp(-EBR/x6)*x5;
        eq(6) = (Ff2*T0+F1*x3-F2*x6)/V2-(kA*exp(-EAR/x6)*x4*dHA+kB*exp(-EBR/x6)*x5*dHB)/Cp+u2/(Cp*V2); 
        eq(7) = (F2*x4-(FD+FR)*(alA*x7/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x7)/V3;
        eq(8) = (F2*x5-(FD+FR)*(alB*x8/(alA*x7+alB*x8+alC*(1-x7-x8)))-F3*x8)/V3;
        eq(9) = (F2*x6-(FD+FR)*x9-F3*x9)/V3+u3/(Cp*V3);

        x1 = x1 + interval*eq(1);
        x2 = x2 + interval*eq(2);
        x3 = x3 + interval*eq(3);
        x4 = x4 + interval*eq(4);
        x5 = x5 + interval*eq(5);
        x6 = x6 + interval*eq(6);
        x7 = x7 + interval*eq(7);
        x8 = x8 + interval*eq(8);
        x9 = x9 + interval*eq(9);
    end
    
       
    x1array(i+1,:) =  x1;
    x2array(i+1,:) =  x2;
    x3array(i+1,:) =  x3;
    x4array(i+1,:) =  x4;
    x5array(i+1,:) =  x5;
    x6array(i+1,:) =  x6;
    x7array(i+1,:) =  x7;
    x8array(i+1,:) =  x8;
    x9array(i+1,:) =  x9;

    t(i+1) = Delta* i;
    delx1=x1-x1s;
    delx2=x2-x2s;
    delx3=x3-x3s;
    delx4=x4-x4s;
    delx5=x5-x5s;
    delx6=x6-x6s;
    delx7=x7-x7s;
    delx8=x8-x8s;
    delx9=x9-x9s;
    
    cost1=cost1+1.5*delx1*delx1+1.5*delx2*delx2+0.002*delx3*delx3+0.00005*u1*u1;
    cost2=cost2+1.2*delx4*delx4+1.2*delx5*delx5+0.0018*delx6*delx6+0.00005*u2*u2;
    cost3=cost3+1.15*delx7*delx7+1.15*delx8*delx8+0.0015*delx9*delx9+0.00005*u3*u3;
    
    % cost1_log = [cost1_log cost1];
    % cost2_log = [cost2_log cost2];
    % cost3_log = [cost3_log cost3];


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

u1array(i,:) = u1pred;
u2array(i,:) = u2pred;
u3array(i,:) = u3pred;
figure(2)
subplot(3,1,1),stairs(t,u1array(:,1),'linewidth',1.5);
ylabel('u_1');
subplot(3,1,2),stairs(t,u2array(:,1),'linewidth',1.5);
ylabel('u_2');
subplot(3,1,3),stairs(t,u3array(:,1),'linewidth',1.5);
ylabel('u_3');
xlabel ('Time(s)')

save ydmpc2
