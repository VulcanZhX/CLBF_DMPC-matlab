function Xns=StatePredict(x0,Uop,dd,A,B,D,P,nd)
% Usage: To Predict the future state of one subsystem
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%   


[m,n]=size(D);

X=x0;
Xns=x0(m-nd+1:m,1);

for i=1:P-1
    X=A*X+B*Uop(i+1)+D*dd(n*i+1:n*(i+1));
    Xns=cat(1,Xns,X(m-nd+1:m,1));
end
