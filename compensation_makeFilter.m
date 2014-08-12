function [B,A]=compensation_makeFilter(w_Sk)
    sampling_freq   = 1/abs(w_Sk(2)-w_Sk(1));
    order           = 2;
    cut_off_freq    = 80;
    peak_to_peak_dB = 0.5;      % This is the Matlab doc recommend first guess
    stopband_atten  = 20;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A]           = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
end