function [Int,t,Ek]=compressor_toTimeDomain(I,l,varargin)
    %Function compressor_toTimeDomain(I,l,varargin)
    %Calculates the Fourier transform for the provided intensity,
    %wavelength and phase.
    %If phase is left out it will be set to constant zero.
    %Optionally a fourth argument can be set. This argument is the factor
    %by which the pahse will be padded with zeros (for higher resolution in
    %the plots).
    %Possible calling methods:
    %   compressor_toTimeDomain(I,l)
    %   compressor_toTimeDomain(I,l,p)
    %   compressor_toTimeDomain(I,l,3)
    %   compressor_toTimeDomain(I,l,p,3)
    if nargin > 4
        error('Too many arguments!\nUsage: compensation_calcFourierlimit intensity wavelength [phase] [padded zeros factor]\n')
    end
    switch nargin
        case 4
            % In this case the caller provided a phase plus a factor for
            % padding zeros
            L = length(l);
            N = L * varargin{2};
            p = [varargin{1}; zeros(N,1)];
            I = [I; zeros(N,1)];
            % We have to take special measures for l because here padding
            % zeros makes no sense. We need to continue the wavelength
            % array
            step = l(end) - l(end-1);
            padd = 1:1:N;
            padd = padd * step + l(end);
            l = [l; padd'];
        case 3
            % In this case the caller provided either a phase or a padding
            % factor. This must be checked and processed accordingly
            if (length(varargin{1}) ~= length(l))
                L = length(l);
                N = L * varargin{1};
                Np = length(l) * (varargin{1}+1);
                I = [I; zeros(N,1)];
                step = l(end) - l(end-1);
                padd = 1:1:N;
                padd = padd * step + l(end);
                l = [l; padd'];
                p = zeros(Np,1);
            else
                p = varargin{1};
            end
        case 2
            % No phase provided -> Fourier limit calculation
            p = 0;
        otherwise
            % Everything that falls into here means a missuse of the
            % function -> issue an error
            error('Too few arguments!\nUsage: compensation_calcFourierlimit intensity wavelength [phase] [padded zeros factor]\n')
    end
    Sk_cplx  = sqrt(I) .* exp(-1i.*p); % For p=0 this equals constant phase -> Fourier limit
    % x = abs(x).*exp(1j.*angle(x))
    [t,Ek]   = iSpeck_Fourier(l.*1e-9,Sk_cplx);
    t        = t.*1e15;
    Int      = abs(trapz(t,abs(Ek).^2));
end