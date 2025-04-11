
clear all;   

% ----- PARAMETERS OF PLANT AND ALL MODELS ---------
%if LinearModel(SamplingTime)==1
   load Plant_Model_Data.mat; 
%end
% ----- OUTPUT REFERENCES -----------
Yr=SetpointCurve(Yr0);

% ----- PARAMETERS OF MPC --------
steprange=130;                                                    % ·ÂÕæ²½³¤                   
P = 10; 

% the weighting matrix
M=P;
q=zeros(P,Ns);
r=zeros(M,Ns);
for i=1:Ns
    q(:,i)=ones(P,1);
    r(:,i)=ones(M,1);
end

% ============= start simulation loop ==============
% data for plant model
TP1=780*ones(ACC.Ny,1);
XX(:,1)=XX0;

% data for Decentralized MPC.
% The states of neighbors

tmp=[];
for i=1:P, tmp=cat(1,tmp,TP1); end 
XN0(:,1)=tmp;

for i=2:Ns
    tmp=[];
    nxi=length(X0{i-1});
    tmp2=X0{i-1}(nxi-ACC.Ny+1:nxi);
    for j=1:P,tmp=cat(1,tmp,tmp2); end
    XN0(:,i)=tmp;
end

% The previous manipulated variable

Uop=ones(Ns,1);
utmp=(Flux/450).^1.3;
for i=1:Ns
    if ACC.header(i)==1
        Uop(i,1)=utmp(i);
    else
        Uop(i,1)=1;
    end
     
end
UUop=ones(P,1)*Uop';


% Start simulation.....
for i=1: steprange
    

% ------------Catalog XN(s)------------------------------------------------
    tmp=[];
    for k=1:P, tmp=cat(1,tmp,TP1); end  
    XN(:,1,i)=tmp;
    for s=2:Ns
        if i==1
            XN(:,s,i)=StatePredict(X0{s-1},UUop(:,s-1),XN0(:,s-1),A{s-1},B{s-1},D{s-1},P,ACC.Ny);
        else
            XN(:,s,i)=StatePredict(X0{s-1},UUop(:,s-1),XN(:,s-1,i-1),A{s-1},B{s-1},D{s-1},P,ACC.Ny);
        end
    end
%-------------Calculate velocity of plate by CONTROLLER-------------------- 
    for s=1:Ns
        if ACC.header(s)==1
           Yrr=Yr(s,i+1:i+P)';
           % Yrr=Yr(s,i+1)*ones(P,1);
           UUop(:,s)=MPCD(A{s},B{s},C{s},D{s},Yr(s,i)*ones(P,1),X0{s},XN(:,s,i),diag(q(:,s)),diag(r(:,s)),P,Uop(s));
        else
           UUop(:,s)=ones(P,1);
        end
        Uop(s,1)=UUop(1,s);
    end
%---------------Calculate the Temperature in Next Step---------------------

    utmp=[];
    for s=1:Ns
        if ACC.header(s)==1
            utmp=cat(1,utmp,Uop(s));
        end        
    end
    ucen(:,i)=utmp; %#ok<AGROW>
   
    XX(:,i+1)=AA* XX(:,i)+BB*ucen(:,i)+DD*[TP1;ones(Ns-sum(ACC.header),1)];
    YY(:,i+1)=CC*XX(:,i+1); %#ok<AGROW>

    ncx=sum(ACC.divnum);
    ncy=ACC.Ny;
    for ix=1:ncx
        TT(:,ix,i)=XX((ix-1)*ncy+1:ix*ncy,i);
    end
    for s=1:Ns
        X0{s}=XX((ACC.subs(s)-1)*ncy+1:ACC.sube(s)*ncy,i);
    end
    UU(:,i)=Uop';
end
% mesh(TT(:,:,i));
save Performance_DMPC.mat