function [J] = PertubedJ(nel,K,theta,nu,th,xx,yy,elem,constraint,Tri,g,cond,force,CD_Level_2)
%%% Function to obtain Jacobian using pertubation
%%% CD_Level_2 = 1; 2nd order finite differencing
%%% CD_Level_2 = 0; 1st order finite differencing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = eps^(1/3); 
nx = nel; % degrees of freedom
nf = size(K,1); 
J = zeros(nf,nx); 
thetan = theta;

for n = 1:nx
    delta = zeros(nx, 1);
    delta(n) = delta(n)+dx;
    x=thetan(n)*ones(nel,1);
    if CD_Level_2
        %%%2nd order central diff%%%%
        dF1 = -FMDL(x+2*delta,nu,th,constraint,Tri,g,cond,force);
        dF2 = 8*FMDL(x+delta,nu,th,constraint,Tri,g,cond,force);
        dF3 = -8*FMDL(x-delta,nu,th,constraint,Tri,g,cond,force);
        dF4 = FMDL(x-2*delta,nu,th,constraint,Tri,g,cond,force);
        dF = dF1+dF2+dF3+dF4;
        J(:, n) = dF(:)/dx/12;
    else
        %%%1st order central fdiff order%%%
        dF = FMDL(x+delta,nu,th,constraint,Tri,g,cond,force)-FMDL(x-delta,nu,th,constraint,Tri,g,cond,force);
        J(:, n) = dF(:)/dx/2; 
    end   
end