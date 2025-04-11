function [A,B,C,D,Tv]=matricesSubsystem(NoSubsys,T_m)
% Usage: To construct the subsystem model of ACC
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%   


global ACC;
global vel;
global dt;
% Initdata;
% T_m  = T(:,sum(ACC.divnum(1:NoSubsys))-ACC.divnum(NoSubsys)+1:sum(ACC.divnum(1:NoSubsys)));
dx   = ACC.dx;
dy   = ACC.dy;
m    = ACC.Ny;
n    = ACC.divnum(NoSubsys);
ctype = ACC.header(NoSubsys);
%function [AM,BV,DM,CV,Tv]=GetSubModel(TM,CType,deltax,deltay,vel,dt,Tinf)
if ctype==1
    Tinf=ACC.Tw;
else
    Tinf=ACC.Tm;
end

Tv=[];
for i=1:n,Tv=cat(1,Tv,T_m(:,i));end
gamma  = vel/dx;            
dim    = m * n;

% Construct AM
% A(i)
A  =  triu(ones(m,m),-1) - triu(ones(m,m),2) - 3 * eye(m,m);
itmp = 2; A(1,2) = itmp; A(m,m-1) = itmp; A(m,m) = -1*itmp; A(1,1) = -1*itmp;

AA=[]; for j=1:n, AA((j-1)*m+1:j*m,(j-1)*m+1:j*m)=A; end

alpha = zeros(dim,1);  for i=1:dim, alpha(i)=fAlpha(Tv(i),dy); end
alpha_m = diag(alpha);

Am = gamma * (triu(ones(dim,dim),-m)-triu(ones(dim,dim),-m+1)- eye(dim,dim));
Am = alpha_m * AA + Am;

% Construct B
B =zeros(dim,1);
for i=1:n
    bata=fBata(Tv((i-1)*m+1),dy);
    B((i-1)*m+1) =-2* bata * transfercoeff(Tv((i-1)*m+1),ctype,Tinf,vel)*(Tv((i-1)*m+1)-Tinf);
    bata=fBata(Tv(i*m)      ,dy);%*dt;
    B(i*m)       =-2* bata * transfercoeff(Tv(i*m)      ,ctype,Tinf,vel)*(Tv(    i*m)-Tinf);
end
% Construct D and C
Dm = zeros(dim,m);
Dm(1:m,1:m)= gamma* eye(m);
C=[zeros(1,dim-m),0.5,ones(1,m-2),0.5]/(m-1); 

% time decrete;
Bm=[B,Dm];

%tic
csys = ss(Am,Bm,C,0);
dsys = c2d(csys,dt);
A=dsys.a;
B=dsys.b(:,1);
C=dsys.c;
D=dsys.b(:,2:m+1);
%toc
% =========================================================================
function alpha=fAlpha(T,dy)
 alpha =tempconductor(T) / dy / dy ;
% 
function bata=fBata(T,dy)
 lumda = conductor(T);
 bata  = fAlpha(T,dy) / lumda * dy;
% 
function h=transfercoeff(T,CType,Tinf,vel)
if CType==0
    h=alpha_aircooling([T(1),T(1)],Tinf);
   % h=h1;
else
    h   =  alpha_water(450,T,Tinf,vel,6);
end
