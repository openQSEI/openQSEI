function [smooth_pr] = WeightedSmoothnessPrior(g,sigvar,corrlength)
%%% Weighted smoothness prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = length(g);
c = (10e-4)*sigvar;
a = sigvar - c;
b = sqrt(-corrlength^2/(2*log(.01)));
smooth_pr = zeros(len,len);

for ii = 1:len
    for jj = ii:len
        dist_ij = norm(g(ii,:)-g(jj,:));
        t_ij = a*exp(-dist_ij^2/(2*b^2));
        if ii == jj
            t_ij = t_ij + c;
        end
        smooth_pr(jj,ii) = t_ij;
        smooth_pr(ii,jj) = t_ij;
    end
end







