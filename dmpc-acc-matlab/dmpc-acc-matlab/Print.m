% Print picture
load Performance_DMPC.mat
%load Performance_MPC.mat
close all
Clinetype={'go-','m--','b+-','r-' ,'k:','k-.'};
mksize=[1.5  4 3 1 1 1]  ;
l=[1 2 3 5];
timeLim=[25 60];
TTmean=[];
Nm=size(TT,1);

for i=1:steprange
    TTmean(i,1) = [0.5 ones(1,Nm-2) 0.5] * TT(:,1,i)/(Nm-1);
    TTmean(i,2:Ns) = [0.5 ones(1,Nm-2) 0.5] * TT(:,4:2:34,i)/(Nm-1);
    TTmean(i,Ns+1) = [0.5 ones(1,Nm-2) 0.5] * TT(:,38,i)    /(Nm-1);
end
x = [1:100];
figure
for s=1:Ns
    plot(x,TTmean(1:100,s)); hold on %#ok<NBRAK>   
    plot(x,Yr(s,1:100)'); hold on     %#ok<NBRAK>
    str=int2str(s);
    ylabel(strcat('Temperature y{\fontsize{6}',int2str(s),'}'));
end


