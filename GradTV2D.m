function grad = GradTV2D(R, theta, beta)
R = sparse(R);
Rx   = R(1:end/2,:);
Ry   = R(end/2+1:2*end/2,:);
Rxtheta = Rx*theta;
Rytheta = Ry*theta;
temp = sqrt((Rxtheta).^2 + (Rytheta).^2  + beta);
grad=R'*[Rxtheta./temp;Rytheta./temp];
