function val=observe_compensation_min_Func(poly,omega,intensity,phase,lambda,fourierlimit)

    %% Too make function definitions easier we define a bunch of global variables

    global I_Sk l_Sk p_Sk w_Sk
    global fit_I_Sk fit_l_Sk fit_p_Sk fit_w_Sk
    global Int_F
    global p order

    %%

    
    %%
    P = polyval([poly 0 0],omega);
    D = phase - P;
    Sk_cplx = sqrt(intensity) .* exp(1i*D);
    [t,E] = Speck_Fourier(lambda',Sk_cplx);
    t = t * 1e15;
    Int = abs(trapz(t,abs(E).^2));
    E = E .* sqrt(fourierlimit ./ Int);
    val = 1-max(abs(E).^2);
end