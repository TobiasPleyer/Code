function val=observe_compensation_min_Func(polynomial)

    %% Too make function definitions easier we define a bunch of global variables

    global I_Sk l_Sk p_Sk w_Sk
    global Int_F

    %%

    
    %%
    P = polyval([polynomial 0 0],w_Sk);
    D = p_Sk - P;
    Sk_cplx = sqrt(I_Sk) .* exp(1i*D);
    [t,E] = Speck_Fourier(l_Sk',Sk_cplx);
    t = t * 1e15;
    Int = abs(trapz(t,abs(E).^2));
    E = E .* sqrt(Int_F ./ Int);
    val = 1-max(abs(E).^2);
end