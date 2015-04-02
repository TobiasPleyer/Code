function [P,GD,GDD]=compressor_GDD_from(wavelength,phase,varargin)
    % Function compressor_GDD_from(wavelength,phase,varargin)
    % This function aims to calculate the GD and GDD values of the provided
    % phase. This is done in the following steps:
    %       1. Transform into frequency domain
    %       2. Fit a polynomial to the phase
    %       3. In frequency domain form the derivatives of that polynomial
    %       4. GD and GDD are the 1st and 2nd of those derivatives
    %Calling options:
    %   compressor_GDD_from(wavelength,phase)
    %   compressor_GDD_from(wavelength,phase,'-')
    %   compressor_GDD_from(wavelength,phase,'+',lower_limit,upper_limit)
    c = 299792458;
    if nargin > 2
        if varargin{1} == '+'
            sign = 1;
        elseif varargin{1} == '-' 
            sign = -1;
        else
            sign = 1;
        end
        if nargin == 5
            lower = varargin{2};
            upper = varargin{3};
        end
    else
        sign = 1;
    end
    % We know the phase in wavelength domain -> transform to frequency space
    % and calculate the GD and GDD
    frequency = (2*pi*c) ./ (wavelength .* 1e-9);
    frequency = fliplr(frequency')'; % The fliplr comes because we change our x-axis: l --> w => ascending --> descending
    phase = fliplr(phase')';
    if nargin == 5
        idx_lower = auxiliary_find_closest_idx(wavelength,lower);
        idx_upper = auxiliary_find_closest_idx(wavelength,upper);
        m = min(idx_lower,idx_upper);
        M = max(idx_lower,idx_upper);
        fit_wavelength = wavelength(m:M);
        fit_phase = phase(m:M);
        fit_frequency = (2*pi*c) ./ (fit_wavelength .* 1e-9);
        fit_frequency = fliplr(fit_frequency')';
        fit_phase = fliplr(fit_phase')';
        poly = polyfit(fit_frequency,sign.*fit_phase,10);
    else
        poly = polyfit(frequency,sign.*phase,10);
    end
    P = polyval(poly,frequency);
    poly_deriv = polyder(poly);
    GD = polyval(poly_deriv,frequency) .* 1e15;
    poly_deriv = polyder(poly_deriv);
    GDD = polyval(poly_deriv,frequency) .* 1e30;
    P = fliplr(P')';
    GD = fliplr(GD')';
    GDD = fliplr(GDD')';
end