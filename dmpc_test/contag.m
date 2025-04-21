clear
a=1060;
b=22;
c=22;
d=0.52;
P=[a b; c d];

F=5;V=1; 
k0=8462639.288942123;
E=50000;
R=8.314;
T0=300; 
Dh=-11500;
rho=1000; 
Cp=0.231; 
%inital values for Q and CA (Assume:steady state)
Qs=0;
CA0s=4;

%steady state: CAs=2      Ts=400
CAs=1.2; Ts=438;

%Lie derivate
%%%simplify the parameters
  A=-F/V;
  C=F*(CA0s-CAs)/V;
  D=Dh/(rho*Cp);
  EE=F*(T0-Ts)/V;


%x1:=CA-CAs         x2:=T-Ts           u1=CA0-CA0s     u2=Q-Qs=Q

%range for x1: 1~3     range for x2:300~500
figure;
axis([0 4 200 600]);
for x1=-1:0.05:4.5
    x1=x1-CAs;
    for x2=250:1:600
        x2=x2-Ts;
LfV=(2*a*x1+(b+c)*x2)*(A*x1-k0*exp(-E/(R*(x2+Ts)))*(x1+CAs)*(x1+CAs)+C)...
    +((b+c)*x1+2*d*x2)*(A*x2-D*k0*exp(-E/(R*(x2+Ts)))*(x1+CAs)*(x1+CAs)+EE);

%LgV=(2*a*x1+(b+c)*x2)*(F/V)+((b+c)*x1+2*d*x2)/(rho*Cp*V);
LgV=((b+c)*x1+2*d*x2)/(rho*Cp*V);


if (LgV ~=0) 
    h2x=-(LfV+sqrt((LfV^2+LgV^4)))/LgV;
else
    h2x=0;
end

%input constraints
if (h2x>5e5) h2x=5e5;
end
if (h2x<-5e5) h2x=-5e5; end



% dV(x)/dt
dVx=(2*a*x1+(b+c)*x2)*(A*x1-k0*exp(-E/(R*(x2+Ts)))*(x1+CAs)*(x1+CAs)+C)...
    +((b+c)*x1+2*d*x2)*(A*x2-D*k0*exp(-E/(R*(x2+Ts)))*(x1+CAs)*(x1+CAs)+EE+h2x/(rho*Cp*V))

if (dVx<0) plot(x1,x2,'b*');hold on; 
%else plot(x1,x2,'r*');
end

    end
end



%Plot ellipse
% hold on;
% syms x y;
% for i=(a-100):20:(a+100)
%     
%   h=ezplot(i*x^2+44*x*y+0.52*y^2-368,[-2 2 -100 100]);
%   set(h,'Color','red')
%   hold on;
% end

