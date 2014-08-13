function [t,E]=compensation_calcCompression(p_orig,p_comp,I_Sk,l_Sk,Int_F)
    d       = p_orig-p_comp;
    [t,E]   = compensation_calcFourierlimit(I_Sk,l_Sk'*1e-9,d);
    t       = t * 1e15;
    Int     = abs(trapz(t,abs(E).^2));
    % Scale the integrals to allow for comparison with Fourier limit
    E       = E .* sqrt(Int_F ./ Int);
end