function Xns=StatePredict2(x0,Uop,dd,A,B,D,P,ny,nx1)
% Usage: To estimate the future state of one subsystem using neighbourhood
% model. the output of this function is the last colum temperature sequence
% of subsystm i.
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%   
%

[m,n]=size(D);

X=x0;
Xns=x0((nx1-1)*ny+1:nx1*ny,1);

for i=1:P-1
    X=A*X+B*Uop(i+1)+D*dd(n*i+1:n*(i+1));
    Xns=cat(1,Xns,X((nx1-1)*ny+1:nx1*ny,1));
end
