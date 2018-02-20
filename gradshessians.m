function [grad_TV,Hess_TV,gradq,Hessq,gradq2,Hessq2,g_smooth,Hess_smooth,cost_params] = gradshessians(theta,cost_params,Fnorm)
%%% Function to obtain the gradients and Hessians
%%% related to the prior models and constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = cost_params.g;
SET_TV_ON = cost_params.SET_TV_ON;
SET_MIN_D = cost_params.SET_MIN_D;
minval      = cost_params.minval;
ii = cost_params.ii;
alpha = cost_params.alpha;
Esim = cost_params.Esim;
%theta = cost_params.theta;
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
tmax = cost_params.tmax;

if SET_TV_ON
    grad_TV = GradTV2D(R,theta,beta);
    Hess_TV = HessTV2D(R,theta,beta);
    grad_TV(isnan(grad_TV))=0;
    Hess_TV(isnan(Hess_TV))=0;
else
    grad_TV = zeros((sNinv),1);
    Hess_TV = zeros((sNinv),(sNinv));
end

if SET_MIN_D
    tminS_1 = minval;
    tminS2_1 = 0.99*minval;
    cost_params.tminS_1 = tminS_1;
    C2 = 1e-4;
    Q = C2*Fnorm(1);
    aS2_1 = -Q/((tminS_1^2-tminS2_1^2)-2*tminS_1*(tminS_1-tminS2_1));
    bS2_1 = -2*aS2_1*tminS_1;
    cS2_1 = -aS2_1*tminS_1^2-bS2_1*tminS_1;
    cost_params.aS2_1 = aS2_1;
    cost_params.bS2_1 = bS2_1;
    cost_params.cS2_1 = cS2_1;
    gradq = zeros((sNinv),1);
    Hessq = zeros((sNinv),(sNinv));
else
    gradq = zeros((sNinv),1);
    Hessq = zeros((sNinv),(sNinv));
end

if SET_MAX_D
    cost_params.tmax = tmax;
    tmax2 = 2*tmax;
    C3 = 1e-4;
    Q = C3*Fnorm(1);
    a2 = -Q/((tmax^2-tmax2^2)-2*tmax*(tmax-tmax2));
    b2 = -2*a2*tmax;
    c2 = -a2*tmax^2-b2*tmax;
    cost_params.a2 = a2;
    cost_params.b2 = b2;
    cost_params.c2 = c2;
    
    gradq2 = zeros((sNinv),1);
    Hessq2 = zeros((sNinv),(sNinv));
else
    gradq2 = zeros((sNinv),1);
    Hessq2 = zeros((sNinv),(sNinv));
end

%%%%%%%%% grads and hessians related to smoothness priors %%%%%%%
%gradients
g_smooth = ((theta - thetaexp(1)*ones(length(thetaexp),1))'*(iGamma))';
%hessians
Hess_smooth = iGamma;