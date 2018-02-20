function sk = linesearch_uniform(theta,dtheta,NMAX,fno,cost_params)
%%% Uniform linesearch to obtain "optimal" stepsize s_k (sk)
%%% for LS optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minval=cost_params.minval;
L_ON=cost_params.L_ON;
kk=1e-10;
Niternonref = 100;
Niterref =100;
S0 =0.02;

if L_ON
    S0 =1;
end
    

%%%% Initial search %%%%
for i=1:Niternonref
    if i==1
    else
        kk = kk+S0;
    end
    T = theta + kk*dtheta;
    ineg = find(T<minval);
    T(ineg)=minval;
    
    F(i)=costfun(T,cost_params);
    
    k(i) = kk;
    
    figure(fno), clf,
    plot(k,F,'b+')
    title('Line-search'), drawnow
    
end

idx = find(F == min(F(:)));

%%%%initial step%%%%
q = k(idx);

%%%%refine%%%%
kkk = q-S0;
if kkk<0
kkk=1e-10;
end

for i=1:Niterref
    kkk = kkk+S0/(Niterref/3);
    T = theta + kkk*dtheta;
    ineg = find(T<minval);
    T(ineg)=minval;
    
    Frefined(i)=costfun(T,cost_params);
    
    kq(i) = kkk;
    
    figure(fno), clf,
    plot(kq,Frefined,'b+')
    title('Line-search'), drawnow
    
end

idx = find(Frefined == min(Frefined(:)));
sk = kq(idx);

figure(fno);
subplot(2,1,1);
semilogy(k,F,'b+'),title('Initial Linesearch'),set(0,'defaulttextInterpreter','latex')
subplot(2,1,2);
semilogy(kq,Frefined,'bo'),title('Refined Linesearch'),set(0,'defaulttextInterpreter','latex')
drawnow