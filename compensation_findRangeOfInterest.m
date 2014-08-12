function [fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk]=compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk)
    M        = max(I_Sk);
    perc     = 0.05;
    fit_w_Sk = w_Sk(I_Sk > perc*M);
    % Now find the fringe values
    lower  = min(fit_w_Sk);
    higher = max(fit_w_Sk);
    % Now we use these values to include everything in between
    fit_w_Sk = w_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_p_Sk = p_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_I_Sk = I_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_l_Sk = l_Sk(w_Sk >= lower & w_Sk <= higher);
end