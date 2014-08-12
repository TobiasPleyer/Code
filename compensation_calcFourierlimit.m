function [Int_F,t_F,Ek_F]=compensation_calcFourierlimit(I_Sk,l_Sk)
    Sk_cplx    = sqrt(I_Sk) .* exp(1i*0); % This equals constant phase -> Fourier limit
    [t_F,Ek_F] = Speck_Fourier(l_Sk*1e-9,Sk_cplx);
    t_F        = t_F * 1e15;
    Int_F      = abs(trapz(t_F,abs(Ek_F).^2));
end