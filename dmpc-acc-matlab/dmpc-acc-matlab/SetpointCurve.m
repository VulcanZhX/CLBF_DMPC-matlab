function [ r ] = SetpointCurve( r0 )
% to initialize the temperature set piont of each subsystem
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university
%   

Ns=size(r0,1);
r=ones(Ns,150);
for i=1:150
    r(:,i)=r0;
end

for k=1:Ns
    if k==1 || k==3 || k==2
        r(k,20:150)=r(k,20:150)+20;
    elseif k==4|| k==12||k==6
        r(k,30:70)=r(k,30:70)+25;
    elseif k==9|| k==5 || k==11
        r(k,20:80)=r(k,20:80)+20;
    else
        r(k,60:150)=r(k,60:150)+20;
    end
end