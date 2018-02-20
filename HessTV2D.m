function H = HessTV2D(R, theta, beta)
R = sparse(R);
Rx   = R(1:end/2,:);
Ry   = R(end/2+1:2*end/2,:);
Rxtht = Rx*theta;
Rytht = Ry*theta;
tmp = ((Rxtht).^2 + (Rytht).^2 + beta);
n=numel(tmp);
 Q=[spdiags(1./sqrt(tmp)-Rxtht.^2./tmp.^1.5,0,n,n) spdiags(-(Rxtht.*Rytht)./tmp.^1.5,0,n,n);
     spdiags(-(Rytht.*Rxtht)./tmp.^1.5,0,n,n) spdiags(1./sqrt(tmp)-Rytht.^2./tmp.^1.5,0,n,n)];   
H=R'*Q*R;