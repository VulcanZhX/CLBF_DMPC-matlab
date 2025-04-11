function Uop=MPCD(A,B,C,D,Yr,x,dd,Q,R,P,upr)
% Usage: Uop=MPCD(A,B,C,D,Yr,x,dd,Q,R,P,upr)
% -------------------------------------------------------------------------
% Function: To slove a MPC optimizaion problem for discreted linear system.
%
% System Model:
%                 x(k+l+1|k) = A * x(k+l|k) + B * u (k+l) + d(k+l|k);
%                 y(k+l+1|k) = C * x(k+l|k);
% Output:
%      Uop = [u(k), u(k+1), ... , u(k+P-1)]'
% Input:
%      dd and Yr is a sequence  eg.   ( [y1(1) y2(1) .......y1(p),y2(p)]')
%      x is initial state of system. and upr is the previous input of
%      system.
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%

G=[]; tmp=C;

m=size(C,1); n=size(B,2);l=size(D,2);
H0=zeros(m,n*P);H=[];
L0=zeros(m,l*P);L=[];

for i=1:P
    H0=[tmp*B,H0(:,1:(P-1)*n)];
    H((i-1)*m+1:i*m,:) = H0;
   
    L0=[tmp*D,L0(:,1:(P-1)*l)];
    L((i-1)*m+1:i*m,:) = L0;
    tmp=tmp*A;
    G=cat(1,G,tmp);
end


% construct coefficient matrix of PQ
Q_bar = H'* Q * H + R; 
f     = -H' * Q * ( Yr - G * x - L * dd );
% Uop=-inv(Q_bar)*f;


% ++++ to construct boundary, constraint coefficient matrix ++++
dimU=n*P;
Ip=eye(dimU);
T=ones(dimU,dimU);
T=eye(dimU) + tril(T,-n-1) - tril(T,-n);
C=cat(1,Ip,-Ip,T,-T);

ulmax=2;  ulmin=0;
dumin=-2;dumax=2;
DUmax=[dumax+upr;dumax*ones(dimU-n,1)];
DUmin=[dumin+upr;dumin*ones(dimU-n,1)];

l   = cat(1,ulmax*ones(dimU,1), -ulmin*ones(dimU,1), DUmax, -DUmin);

% using quadratic programing to optimal constraint MPC problem
Uop = quadprog(Q_bar,f,C,l);
% % %--------------------------------------------------------------------------
