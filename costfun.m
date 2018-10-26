function [Fnew,usim] = costfun(E_GN,cost_params)
%%% Cost function (l) updated in each iteration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% To avoid confusion of variables "theta" is called "E_GN"

g = cost_params.g;
SET_TV_ON = cost_params.SET_TV_ON;
minval      = cost_params.minval;
alpha = cost_params.alpha;
um = cost_params.um;
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
thetaexp=cost_params.thetaexp;
force = cost_params.force;
SET_MAX_D = cost_params.SET_MAX_D;
SET_MIN_D = cost_params.SET_MIN_D;
tmax = cost_params.tmax;
lambda = cost_params.lambda;
Ln = cost_params.Ln;

minind = find(E_GN < minval);
q = sum(aS2_1*(E_GN(minind)).^2 + bS2_1*E_GN(minind) + cS2_1);

if SET_MAX_D
    a2 = cost_params.a2 ;
    b2 = cost_params.b2;
    c2 = cost_params.c2;
    maxind = find(E_GN > tmax);
    q2 = sum(a2*(E_GN(maxind)).^2 + b2*E_GN(maxind) + c2);
else
    q2=0;
end

[usim,K] = FMDL(E_GN,nu,th,constraint,Tri,g,cond,force);

if SET_TV_ON  
    Rt = R*E_GN;
    nanRt = isnan(Rt);
    Rt(nanRt)=0;
    Rx = Rt(1:end/2,:);
    Ry = Rt(1+end/2:end,:);
    t = sum(Ai.*sqrt( Rx.^2 + Ry.^2 + beta ));
    Fnew = 0.5*norm(Ln*(um-usim))^2  + q + q2 + alpha*t;
end

if L_ON
    Fnew = 0.5*norm(Ln*(um-usim))^2 + 0.5*lambda*norm(E_GN) + q + q2;
end

if FIRST_ORDER
    smooth = 0.5*(E_GN - thetaexp(1)*ones(numel(E_GN),1))'*iGamma*(E_GN - thetaexp(1)*ones(numel(E_GN),1));
    Fnew = 0.5*norm(Ln*(um-usim))^2  + q + q2 + smooth;
end
