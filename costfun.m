function [Fnew,usim] = costfun(sig,cost_params)
%%% Cost function (l) updated in each iteration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = cost_params.g;
SET_TV_ON = cost_params.SET_TV_ON;
minval      = cost_params.minval;
ii = cost_params.ii;
alpha = cost_params.alpha;
Esim = cost_params.Esim;
theta = cost_params.theta;
um = cost_params.um;
sNinv = cost_params.sNinv;
thetaexp = cost_params.thetaexp;
a = cost_params.a;
aS2_1 = cost_params.aS2_1 ;
bS2_1 = cost_params.bS2_1;
cS2_1 = cost_params.cS2_1;
R = cost_params.R;
Ai = cost_params.Ai;
beta = cost_params.beta;
nu = cost_params.nu;
th = cost_params.th;
iGamma = cost_params.iGamma;
FIRST_ORDER = cost_params.FIRST_ORDER;
L_ON = cost_params.L_ON;
constraint = cost_params.constraint;
Tri = cost_params.Tri;
cond = cost_params.cond;
force = cost_params.force;
SET_MAX_D = cost_params.SET_MAX_D;
SET_MIN_D = cost_params.SET_MIN_D;
tmax = cost_params.tmax;
lambda = cost_params.lambda;

minind = find(sig < minval);
q = sum(aS2_1*(sig(minind)).^2 + bS2_1*sig(minind) + cS2_1);

if SET_MAX_D
    a2 = cost_params.a2 ;
    b2 = cost_params.b2;
    c2 = cost_params.c2;
    maxind = find(sig > tmax);
    q2 = sum(a2*(sig(maxind)).^2 + b2*sig(maxind) + c2);
else
    q2=0;
end

[usim,K] = FMDL(sig,nu,th,constraint,Tri,g,cond,force);

if SET_TV_ON  
    Rt = R*sig;
    nanRt = isnan(Rt);
    Rt(nanRt)=0;
    Rx = Rt(1:end/2,:);
    Ry = Rt(1+end/2:end,:);
    t = sum(Ai.*sqrt( Rx.^2 + Ry.^2 + beta ));
    Fnew = 0.5*norm((um-usim))^2  + q + q2 + alpha*t;
end


if L_ON
    Fnew = 0.5*norm((um-usim))^2 + 0.5*lambda*norm(sig);
end

if FIRST_ORDER
    smooth = 0.5*(sig - thetaexp(1)*ones(numel(sig),1))'*iGamma*(sig - thetaexp(1)*ones(numel(sig),1));
    Fnew = 0.5*norm((um-usim))^2  + q + + q2 + smooth;
end




