function [EhomF] = BestHomogeneousE(Ehomguess,nu,th,xx,yy,elem,Uel2,Estep,constraint,Tri,g,cond,force)
%%% Function: Best Homogeneous Estimate for E
%%% One parameter LS minimization using a simple confidence
%%% interval search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ehom=[];
ut=[];
Ehom(:,1) = Ehomguess;
idx=[];
for i=2:200
    
    Ehom(:,i) = Ehom(:,i-1) + Estep;
    [u,K] = FMDL(Ehom(:,i),nu,th,constraint,Tri,g,cond,force);
    N1(i)=norm(u-Uel2);
    ut(:,i) = u; 
    idx(i) = i;
    
end
N1 = N1(10:end);
[m,idxx] = min(N1);

EhomF = Ehom(:,idxx);
figure(13)
plot(N1),title('Best Homogeneous Guess'),set(0,'defaulttextInterpreter','latex');

