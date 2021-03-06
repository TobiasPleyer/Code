function val=compensation_minFuncForBruteForce(polynomial)

    %% Too make function definitions easier we define a bunch of global variables

    global I_Sk p_Sk w_Sk
    global Int_F
 
    %%
    P          = polyval(polynomial,w_Sk);
    % In order to increase the sampling rate in the fourier domain, it is
    % necessary to add more information in the spectral domain.
    % Since we cannot invent a phase in this spectral region, a common
    % trick is to add zeros.
    D          = p_Sk - P;
    L          = length(D);
    [w,I,p,l]  = compensation_extendPhaseByZeros(w_Sk,I_Sk,D,L);
    Sk_cplx    = sqrt(I) .* exp(1i*p);
    [t,E]      = Speck_Fourier(l,Sk_cplx);
    t          = t * 1e15;
    Int        = abs(trapz(t,abs(E).^2));
    E          = E .* sqrt(Int_F ./ Int);
    val        = 1-max(abs(E).^2);
end