function [uxsim,uysim] = uxuy(usim)   
%%% Function to obtain x and y displacement components
%%% ux and uy, respectively 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.2.2018 Danny Smyl
%%% Aalto University, Espoo, Finland
%%% Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

countx=0;
county=0;
    for gg = 1:length(usim)
        if mod(gg,2)==1
            countx = countx +1;
            uysim(countx)=usim(gg);
        else
            county = county +1;
            uxsim(county)=usim(gg);
        end
    end