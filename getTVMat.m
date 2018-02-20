function [R,Ai] = getTVMat(g,H,nel)
%%% Function: Get TV parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn=size(g,1);
ne=size(H,1);
R = sparse(2*ne, nn);
Ai = zeros(ne,1);
M = [ -1 0 1;-1 1 0];

for ii=1:ne
    
    gg =g(H(ii,:),:);
    Ao = ((gg(3,1)-gg(1,1)) * (gg(2,2)-gg(1,2))) - (((gg(3,2)-gg(1,2)) * (gg(2,1)-gg(1,1))));
    
    tr = 1/Ao*[gg(2,2)-gg(1,2), gg(1,2)-gg(3,2);...
        gg(1,1)-gg(2,1) gg(3,1)-gg(1,1)];
    
    a = sqrt(sum((gg(2,:)-gg(1,:)).^2));
    b = sqrt(sum((gg(3,:)-gg(2,:)).^2));
    c = sqrt(sum((gg(3,:)-gg(1,:)).^2));   
    s = (a + b + c)/2;
    Ai(ii) = sqrt(s*(s-a)*(s-b)*(s-c));

    R([ii end/2+ii],H(ii,:))=tr*M;   
end

%%% No weighting with respect to element sizes %%%
Ai = ones(size(Ai));



