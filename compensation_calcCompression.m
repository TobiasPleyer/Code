function [t,E]=compensation_calcCompression(p_orig,p_comp,I_Sk,l_Sk,Int_F)
    d       = p_orig-p_comp;
    Sk_cplx = sqrt(I_Sk) .* exp(1i*d);
    [t,E]   = Speck_Fourier(l_Sk'*1e-9,Sk_cplx);
    t       = t * 1e15;
    Int     = abs(trapz(t,abs(E).^2));
    % Scale the integrals to allow for comparison with Fourier limit
    E       = E .* sqrt(Int_F ./ Int);
end