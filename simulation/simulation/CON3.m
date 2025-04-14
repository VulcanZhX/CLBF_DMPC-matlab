function [u3exit,u3,fval9] = CON3(u3s,x7,x8,x9,x7s,x8s,x9s,x4,x5,x6,x4s,x5s,x6s,Delta,interval,pred,u3pred,F2,F3,FD,FR,V3,alA,alB,alC,Cp)
%% Objective Function
fun = @(u3) SUB3(u3,u3s,x7,x8,x9,x7s,x8s,x9s,x4,x5,x6,Delta,interval,pred,F2,F3,FD,FR,V3,alA,alB,alC,Cp);

%% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon =@(u3) nonlinear3(u3,x4,x5,x6,x7,x8,x9,x4s,x5s,x6s,x7s,x8s,x9s,F2,F3,FD,FR,Cp,V3);
%% Bounds (lb <= u <= ub)

lb = [-12.6e5-0.55e5];
ub = [12.6e5+0.55e5];

%% Initial Guess
x0 = u3pred;

%% solve
problem = createOptimProblem('fmincon','objective',fun,'lb',lb,'ub',ub,'x0',x0,'nonlcon',nlcon);   
[u3,fval9,exitflag,info] = fmincon(problem);

%% fval
u3exit(1)= u3(1);

