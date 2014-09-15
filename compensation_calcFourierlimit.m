function [Int,t,Ek]=compensation_calcFourierlimit(I,l,varargin)
    if nargin > 3
        error('Too many arguments!\nUsage: compensation_calcFourierlimit intensity wavelength [phase]\n')
    end
    if nargin == 3
        p = varargin{1};
    else
        p = 0;
    end
    Sk_cplx  = sqrt(I) .* exp(-1i.*p); % For p=0 this equals constant phase -> Fourier limit
    [t,Ek]   = Speck_Fourier(l.*1e-9,Sk_cplx);
    t        = t.*1e15;
    Int      = abs(trapz(t,abs(Ek).^2));
end