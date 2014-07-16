function t_rms = RMS_moritz(t, I)
    I_n = I/(trapz(I));
    
    t_rms = sqrt(trapz(t.^2.*I_n) - trapz(t.*I_n));
end