function [w,GDD]=auxiliary_GDD_from_with_Filter(wavelength,phase,cut_off)
    
    c = 299792458;
    frequency = (2*pi*c) ./ (wavelength .* 1e-9);
    frequency = fliplr(frequency')'; % The fliplr comes because we change our x-axis: l --> w => ascending --> descending
    phase = fliplr(phase')';
    
    [B,A] = makeFilter(frequency,cut_off);
    [w,D] = diffdiff(frequency,phase,B,A);
    [w,D] = diffdiff(w,D,B,A);
    w = 1e9 .* (2*pi*c) ./ fliplr(w')';
    GDD = fliplr(D')' .* 1e30;
end

function [B,A]=makeFilter(x,cut_off)
    sampling_freq   = 1/abs(x(2)-x(1)); % this is slightly incorrect because the frequency is not equi-distant
    order           = 2;
    cut_off_freq    = 1e-9;
    peak_to_peak_dB = 0.5;      % This is the Matlab doc recommend first guess
    stopband_atten  = 20;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
    normed_cutoff   = cut_off;
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A]           = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
end

function [x,D]=diffdiff(x,y,B,A)
    D = diff(y) ./ (x(1)-x(2));
    % Immediately apply our lowpass filter
    D = filtfilt(B,A,D);
    x = x(1:end-1);
end