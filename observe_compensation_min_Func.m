function val=observe_compensation_min_Func(poly,omega,intensity,phase,lambda,fourierlimit)
    P = polyval([poly 0 0],omega);
    D = phase - P;
    Sk_cplx = sqrt(intensity) .* exp(1i*D);
    [t,E] = Speck_Fourier(lambda',Sk_cplx);
    t = t * 1e15;
    Int = abs(trapz(t,abs(E).^2));
    E = E .* sqrt(fourierlimit ./ Int);
    val = 1-max(abs(E).^2);
end