function Y=LinearModel()
% Usage: to calcualte the global model of entire system and the model of
% each susystem, as well as the model of each neighbourhood systems.
% -------------------------------------------------------------------------
% Explaination:
%    Y is logical number to judge if the function is running success 
%    The file of "Plant_Model_Data.mat" inludes the data of Plant in
%    equilibria and the coefficents of models, as well as the initial state
%    and output of each model. 
%
%                                       Yi ZHENG  ¡ò2009.08.17 Version: 1.0 
%                                             Shanghai Jiao Tong university  
%     


% ---------- CALCULATE PLANT INITIAL DATA -------------

%if EquilibriaCalculation(SamplingTime)==1
if EquilibriaCalculation==1    
   load Plant_Data.mat; 
end

Y=0;

% ------------- MODELS INITIALIZATION ----------------
% 1) Model of entire system
  [AA, BB, CC,DD,XX0]=matricesSystem(ACC.T);
  
% 2) Model of each subsystem
for i=1:Ns  
    [A{i}, B{i}, C{i},D{i},X0{i}]=matricesSubsystem(i,ACC.T(:,ACC.subs(i):ACC.sube(i))); 
end

% 3) Model of output-neighborhood
for i=1:Ns-1
    [An{i},Bn{i},Cn{i},Dn{i},Xn0{i}]=matricesNeighborhood(A{i}, B{i}, C{i},D{i},A{i+1},B{i+1},C{i+1},D{i+1},X0{i},X0{i+1});
end

An{Ns}=A{Ns};Bn{Ns}=B{Ns};Cn{Ns}=C{Ns};Dn{Ns}=D{Ns};Xn0{Ns}=X0{Ns};

% equilibria temperature profile along the cooling section 

Yr0 = [0.5 ones(1,ACC.Ny-2) 0.5]*ACC.T(:,ACC.sube)/(ACC.Ny-1);
Yr0 = Yr0';                                                                % initial average temperature 
% Save data to workspace file
clear i;
save Plant_Model_Data.mat;
Y=1;