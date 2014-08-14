function [t,E]=compensation_calcCompression(p_orig,p_comp,I_Sk,l_Sk,Int_F)
    d       = p_orig-p_comp;
    [Int,t,E] = compensation_calcFourierlimit(I_Sk,l_Sk,d);
    % Scale the integrals to allow for comparison with Fourier limit
    E       = E .* sqrt(Int_F ./ Int);
end