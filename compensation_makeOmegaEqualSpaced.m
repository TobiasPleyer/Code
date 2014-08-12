function [l0,w0,w_Sk,I_Sk,p_Sk,l_Sk]=compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk)
    c    = 299792458;
    l0   = 1.030;
    w0   = 2*pi*0.299792458 / l0;
    w_Sk = 2*pi*c ./ (l_Sk*1e-9);
    w_Sk = w_Sk * 1e-15;
    
    % For differentiation it is good to have the x-values equally spaced
    w_Sk2 = linspace(w_Sk(1),w_Sk(end),length(w_Sk));
    I_Sk = interp1(w_Sk,I_Sk,w_Sk2);
    p_Sk = interp1(w_Sk,p_Sk,w_Sk2);
    w_Sk = w_Sk2; clear w_Sk2
    l_Sk = 2*pi*c ./ (w_Sk*1e15);
    l_Sk = l_Sk' * 1e9;
end