% return FT of spectrum
function [t1,E] = Speck_Fourier(lambda, S)
    
    c = 299792458;
    N = length(lambda);
    f = c./lambda;
    Sw = lambda.^2./(2*pi*c).*S; % This possibly leads to problems if this factor has been multiplied beforehand --> check
    
    Fi = linspace(f(1), f(end), N);
    Sw = interp1(f, Sw, Fi, 'cubic');
    
    T = 1/(Fi(1)-Fi(2));
    t1 = -T/2.*linspace(-1,1,N);

    Sw = Sw./max(Sw);
    
    %E = fftshift(fft(ifftshift(Sw))); OLD LINE
    E = fftshift(ifft(ifftshift(Sw)));
    E = E./max(E);

end