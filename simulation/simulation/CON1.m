function [u1exit,u1,fval3] = CON1(u1s,x1,x2,x3,x1s,x2s,x3s,x7,x8,x9,x7s,x8s,x9s,Delta,interval,pred,u1pred,F1,FR,Ff1,V1,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0)
%% Objective Function
fun = @(u1) SUB1(u1,u1s,x1,x2,x3,x1s,x2s,x3s,x7,x8,x9,Delta,interval,pred,F1,FR,Ff1,V1,alA,alB,alC,kA,kB,EAR,EBR,dHA,dHB,Cp,T0,xA0);

%% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon =@(u1) nonlinear1(u1,x1,x2,x3,x7,x8,x9,x1s,x2s,x3s,x7s,x8s,x9s,F1,FR,Ff1,V1,kA,kB,EAR,EBR,dHA,dHB,Cp,T0);
%% Bounds (lb <= u <= ub)

lb = [12.6e5-0.55e6];
ub = [12.6e5+0.55e6];

%% Initial Guess
x0 = u1pred;

%% solve
problem = createOptimProblem('fmincon','objective',fun,'lb',lb,'ub',ub,'x0',x0,'nonlcon',nlcon);   
[u1,fval3,exitflag,info] = fmincon(problem);

%% fval
u1exit(1)= u1(1);

