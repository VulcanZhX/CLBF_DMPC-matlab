function [u2exit,u2,fval6] = CON2(u2s,x4,x5,x6,x4s,x5s,x6s,x1,x2,x3,x1s,x2s,x3s,Delta,interval,pred,u2pred,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
%% Objective Function
fun = @(u2) SUB2(u2,u2s,x4,x5,x6,x4s,x5s,x6s,x1,x2,x3,Delta,interval,pred,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);

%% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon =@(u2) nonlinear2(u2,x1,x2,x3,x4,x5,x6,x1s,x2s,x3s,x4s,x5s,x6s,F1,F2,Ff2,V2,kA,kB,EAR,EBR,dHA,dHB,Cp,T0);
%% Bounds (lb <= u <= ub)

lb = [16.2e5-0.55e5];
ub = [16.2e5+0.55e5];

%% Initial Guess
x0 = u2pred;

%% solve
problem = createOptimProblem('fmincon','objective',fun,'lb',lb,'ub',ub,'x0',x0,'nonlcon',nlcon);   
[u2,fval6,exitflag,info] = fmincon(problem);

%% fval
u2exit(1)= u2(1);

