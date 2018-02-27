function grad = GradTV2D(R, theta, beta)
R = sparse(R);
Rx   = R(1:end/2,:);
Ry   = R(end/2+1:2*end/2,:);
Rxtht = Rx*theta;
Rytht = Ry*theta;
tmp = sqrt((Rxtht).^2 + (Rytht).^2  + beta);
grad=R'*[Rxtht./tmp;Rytht./tmp];
