function [uexit,u,fval] = controller(us,x,xs,Delta,interval,pred,upred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
%% Objective Function
fun = @(u) SUB(u,us,x,xs,Delta,interval,pred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);

%% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon =@(u) nonlinear(u,us,x,xs,Delta,interval,pred,F1,F2,F3,FD,FR,Ff1,Ff2,V1,V2,V3,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);
%% Bounds (lb <= u <= ub)

lb = [(12.6e5-0.55e5)*ones(pred,1);(16.2e5-0.55e5)*ones(pred,1);(12.6e5-0.55e5)*ones(pred,1)];
ub = [(12.6e5+0.55e5)*ones(pred,1);(16.2e5+0.55e5)*ones(pred,1);(12.6e5+0.55e5)*ones(pred,1)];

%% Initial Guess
x0 = upred;

%% solve
problem = createOptimProblem('fmincon','objective',fun,'lb',lb,'ub',ub,'x0',x0,'nonlcon',nlcon);   
[u,fval,exitflag,info] = fmincon(problem);

%% fval
uexit(1)= u(1);
uexit(2)=u(11);
uexit(3)=u(21);

