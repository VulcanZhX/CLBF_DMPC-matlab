function [A,B,C,D,Tv]=matricesSystem(T)
% To construct the overall model of ACC
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%   
%Initdata;
global ACC;
global vel;
global dt;

%T=ACC.T;


dx   = ACC.dx;
dy   = ACC.dy;
m    = ACC.Ny;
n    = sum(ACC.divnum);

ctype=[];
for i=1:length(ACC.divnum);
    ctype(1,ACC.subs(i):ACC.sube(i))=ACC.header(i);
end

% ++ Tv ++
Tv=[];
for i=1:n,Tv=cat(1,Tv,T(:,i));end

% ++ Construct AM ++
gamma  = vel/dx;            
dim    = m * n;

% A
A  =  triu(ones(m,m),-1) - triu(ones(m,m),2) - 3 * eye(m,m);
itmp = 2; A(1,2) = itmp; A(m,m-1) = itmp; A(m,m) = -1*itmp; A(1,1) = -1*itmp;

% AA
AA=[]; for j=1:n, AA((j-1)*m+1:j*m,(j-1)*m+1:j*m)=A; end

% alpha matrice
alpha = zeros(dim,1);  
for i=1:dim 
    alpha(i)=fAlpha(Tv(i),dy); 
end

alpha_m = diag(alpha);

Am = gamma * (triu(ones(dim,dim),-m)-triu(ones(dim,dim),-m+1)- eye(dim,dim));
Am = alpha_m * AA + Am;

% Construct B

nu = sum(ACC.header);
Ns = length(ACC.header);
Bm = zeros(dim,nu);
Dm = zeros(dim,Ns-nu+m);
k=0;
l=0;
for i=1:Ns
    if ACC.header(i)==1
        k=k+1;
        for j=ACC.subs(i):ACC.sube(i)
            Tinf=ACC.Tw;
            bata=fBata(Tv((j-1)*m+1),dy);
            Bm((j-1)*m+1,k)=-2* bata * transfercoeff(Tv((j-1)*m+1),1,Tinf,vel)*(Tv((j-1)*m+1)-Tinf);
            bata=fBata(Tv(j*m)      ,dy);      %  *dt;
            Bm(j*m      ,k)=-2* bata * transfercoeff(Tv(j*m)      ,1,Tinf,vel)*(Tv(    j*m)-Tinf);
        end
    else
        l=l+1;
        for j=ACC.subs(i):ACC.sube(i)
            Tinf=ACC.Tm;
            bata=fBata(Tv((j-1)*m+1),dy);
            Dm((j-1)*m+1,l+m)=-2* bata * transfercoeff(Tv((j-1)*m+1),0,Tinf,vel)*(Tv((j-1)*m+1)-Tinf);
            bata=fBata(Tv(j*m)      ,dy);      %  *dt;
            Dm(j*m      ,l+m)=-2* bata * transfercoeff(Tv(j*m)      ,0,Tinf,vel)*(Tv(    j*m)-Tinf);
        end
    end
end
% Construct D and C
% Dm = zeros(dim,m);
Dm(1:m,1:m)= gamma* eye(m);
% 
% C=[zeros(1,dim-m),1,zeros(1,m-2),0;
%    zeros(1,dim-m),0,zeros(1,m-2),1]; 
%C=[];
%ctemp1=[0.5, zeros(1,m-2),0.5]/(m-1);
C=zeros(Ns,dim);
for s=1:Ns
    C(s,m*(ACC.sube(s)-1)+1: m*(ACC.sube(s)))=[0.5, ones(1,m-2),0.5]/(m-1);
   
%    C=cat(1,C,cat(2, zeros(1,m*(ACC.sube(s)-1)),ctemp1,zeros(1,m*(ACC.sube(s)-1))   ));
end
% time decrete;
BB=[Bm,Dm];
csys = ss(Am,BB,C,0);
dsys = c2d(csys,dt);
A=dsys.a;
B=dsys.b(:,1:nu);
C=dsys.c;
D=dsys.b(:,nu+1:Ns+m);

%figure;
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
else
    h   =  alpha_water(450,T,Tinf,vel,6);
end


