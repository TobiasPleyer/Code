function [new_phase,grating_phase,GDD,TOD] = compensation_applyGrating(omega,phase,varargin)
    % The formulas used below are taken from the script of Prof. Dr. Karsch,
    % pages 31-34 and apply to a grating compressor.

    strUsage = 'Usage:\napplyGrating(omega,phase)\nor\napplyGrating(omega,phase,G,AOI,N)\n';
    strExplanation = 'Where\nomega: Angular frequency vector\nphase: Phase in angular frequency domain\nG: Distance of the gratings [m] (perpendicular)\nAOI: Angle of incidence [rad]\nN: Line number of the grating (d^-1) [#/mm]\n';
    
    if nargin < 2
        fprintf(strUsage)
        fprintf(strExplanation)
    elseif nargin == 2
        G = 0.15;
        AOI = 2*pi * 45/360;
        N  = 300 * 1e3;
        d = 1/N;
    elseif nargin == 5
        G = varargin{1};
        AOI = varargin{2};
        N = varargin{3} * 1e3;
        d = 1/N;
    else
        fprintf(strUsage)
        fprintf(strExplanation)
    end
    
    c = 299792458;
    lambda = 2*pi*c ./ omega;
    
    beta = asin(N.*lambda ./ sin(AOI));
    theta = AOI - beta;
    
    grating_phase = 2*G .* (omega./c.*((1+cos(theta)./cos(beta))) - tan(beta).*(2*pi)./d);
    
    new_phase = phase - grating_phase;
    
    lambda0 = 1.03e-6;
    beta0 = asin(N*lambda0-sin(AOI));
    D1 = G./cos(beta0);
    D2 = -2*D1./c.*lambda0./(2*pi*c) .* (lambda0.*N./cos(beta0)).^2;
    D3 = -3.*D2.*lambda0./(2*pi*c) .* (1 + lambda0.*N.*sin(beta0)./(cos(beta0)).^2);
    GDD = D2 / 1e-30;
    TOD = D3 / 1e-45;
end