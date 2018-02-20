function [g,x,y,Ex,Ey,Tri,constraint,nel,nelx,nely] = meshgen(Lx,Ly,elwidth)
%%% Function to gererate meshing parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%scalable constraint related to the BC's of the FEM%%%
constraint =1e5;

xx =0:elwidth:Lx;
yy =0:elwidth:Ly;
[X,Y] = meshgrid(xx,yy);
g = [X(:),Y(:)];
ginv=g;
x=g(:,1);
y=g(:,2);
Tri = delaunayTriangulation(g(:,1),g(:,2)); %Tri is H
Eloc = incenter(Tri);
Ex= Eloc(:,1);
Ey=Eloc(:,2);
nel=length(Ex);
nelx= Lx/elwidth;
nely= Ly/elwidth;