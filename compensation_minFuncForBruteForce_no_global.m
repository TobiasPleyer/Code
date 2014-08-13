function val=compensation_minFuncForBruteForce_no_global(polynomial,w,I,p,Int_F)

    P          = polyval(polynomial,w);
    % In order to increase the sampling rate in the fourier domain, it is
    % necessary to add more information in the spectral domain.
    % Since we cannot invent a phase in this spectral region, a common
    % trick is to add zeros.
    D          = p - P;
    disp(min(abs(D(:))))
    L          = length(D);
    [w,I,D,l]  = compensation_extendPhaseByZeros(w,I,D,L);
    Sk_cplx    = sqrt(I) .* exp(1i*D);
    [t,E]      = Speck_Fourier(l.*1e-9,Sk_cplx);
    t          = t * 1e15;
    Int        = abs(trapz(t,abs(E).^2));
    E          = E .* sqrt(Int_F ./ Int);
    val        = 1-max(abs(E).^2);
end