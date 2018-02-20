function [u,K] = FMDL(E,nu,th,constraint,Tri,g,cond,force)
%%% Piecewise linear plane stress FEM Forward Model (FMDL)
%%% to obtain displacemnt field u
%%% E,u,nu,th = inhomogeneous elastic modulus, displacement field, homgenous Poisson ratio, plate thickness
%%% Tri, g = triangularization, node points
%%% constraint, force = nodal constraints, force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=g(:,1);
y=g(:,2);

elem = Tri.ConnectivityList;
nel=length(elem(:,1));
np=2*length(x);

Eloc = incenter(Tri);
Ex= Eloc(:,1);
Ey=Eloc(:,2);

elem =[elem, ones(nel,1), ones(nel,1), zeros(nel,1)];
K=sparse(np,np);
F=zeros(np,1);

for i=1:nel
    nod1=elem(i,1);
    nod2=elem(i,2);
    nod3=elem(i,3);
    
    tm  =elem(i,4);
    tp  =elem(i,5);
    
    id=[2*nod1-1  2*nod1  2*nod2-1  2*nod2  2*nod3-1  2*nod3];
    
    E1=E(i)/(1-nu*nu);
    G=E(i)/2/(1+nu);
    C=[E1  nu*E1   0
        nu*E1     E1   0
        0      0    G ];
    
    y32=y(nod3)-y(nod2);
    y13=y(nod1)-y(nod3);
    y21=y(nod2)-y(nod1);
    
    x23=x(nod2)-x(nod3);
    x31=x(nod3)-x(nod1);
    x12=x(nod1)-x(nod2);
    
    Ael=abs((x12*y13-x31*y21))/2;
    
    B=[ y32     0   y13     0   y21     0
        0   x23     0   x31     0   x12
        x23   y32   x31   y13   x12   y21 ]/2/Ael;
    
    DB=C*B;
    kelem=Ael*th*B'*DB;
    K(id,id)=K(id,id)+kelem;
end

for i=1:length(force(:,1))
    z1 =force(i,1);
    dir=force(i,2);
    f =force(i,3);
    l =2*(z1-1)+dir;
    F(l)=F(l)+f;
end

q=max(diag(K))*constraint;

for i=1:length(cond(:,1))
    z2 = cond(i,1);
    dir = cond(i,2);
    l=2*(z2-1) + dir;
    K(l,l)=q;
end

u=K\F;
