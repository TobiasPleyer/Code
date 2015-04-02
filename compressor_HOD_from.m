function [P,HOD]=compressor_HOD_from(wavelength,phase,order,fit_order)
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
    % We know the phase in wavelength domain -> transform to frequency space
    % and calculate the GD and GDD
    frequency = (2*pi*c) ./ (wavelength .* 1e-9);
    frequency = fliplr(frequency')'; % The fliplr comes because we change our x-axis: l --> w => ascending --> descending
    phase = fliplr(phase')';
    poly = polyfit(frequency,phase,fit_order);
    P = polyval(poly,frequency);
    for i=1:order
        poly = polyder(poly);
    end
    HOD = polyval(poly,frequency) .* 10^(order*15);
    P = fliplr(P')';
    HOD = fliplr(HOD')';
end