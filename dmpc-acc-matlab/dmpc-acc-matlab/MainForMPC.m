
clear all;   

% ----- PARAMETERS OF PLANT AND ALL MODELS -------
%if LinearModel==1
   load Plant_Model_Data.mat; 
% end
Yr=SetpointCurve(Yr0);

% ----- PARAMETERS FOR MPC -----
P = 10; 
q=ones(P,Ns);
r=ones(P,sum(ACC.header));
qq=[];rr=[];
for i=1:P
    qq=cat(2,qq,q(i,:));
    rr=cat(2,rr,r(i,:));
end
Q=diag(qq);
R=diag(rr);

% ----- PARAMETERS FOR SIMULATION -----

% initial state and disturbance of entire system

TP1=780*ones(ACC.Ny,1);
XX(:,1)=XX0;

% initial input of system
upre=[];
utmp=(Flux/450).^1.3;
for i=1:Ns
    if ACC.header(i)==1
        upre=cat(1,upre,utmp(i));
    end
end
nu=length(upre);

% ============= start simulation loop ==================
steprange=130;          % ·ÂÕæ²½³¤  
for i=1: steprange
    % calcalte the distrubance of system
    dair(:,1)=ones(Ns-sum(ACC.header),1);
    dall=[TP1;dair];
    dd=[];
    for k=1:P 
        dd=cat(1,dd,dall);
    end
    
    % calculate the output references
    Yrr=[];

    for k=1:P 
        %Yrr=cat(1,Yrr,Yr(:,i+1));   
        Yrr=cat(1,Yrr,Yr(:,i+k));
    end
    %----------Calculate velocity of plate by CONTROLLER--------------- 
    % y1(1) y2(1) .......y1(p),y2(p)
    uuop=MPCD(AA,BB,CC,DD,Yrr,XX(:,i),dd,Q,R,P,upre);
    uop=uuop(1:nu,1);
    upre=uop;
   
    
    
   %-----------Calculate the Temperature in Next Step-------------
   % Flux=Flux+3*rand(1,Ns)-1.5;  % For experiment
    ncx = sum(ACC.divnum);
    ncy = ACC.Ny;
    for ix = 1:ncx
        TT(:,ix,i) = XX((ix-1)*ncy+1:ix*ncy,i);
    end
   
    XX(:,i+1)=AA* XX(:,i)+BB*uop+DD*dall;
    YY(:,i+1)=CC*XX(:,i+1); %#ok<AGROW>
    
    % record input variables
    k=0;
    Uop=ones(1,Ns);
    for s=1:Ns
        if ACC.header(s)==1
            k=k+1;
            Uop(s)=uop(k);
        end
    end
    UU(:,i)=Uop';
 
end
save Performance_MPC.mat
