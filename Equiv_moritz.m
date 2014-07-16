function t_equiv = Equiv_moritz(t, I)
    I_max = max(I);
    
    t_equiv = 1/I_max .* trapz(I);
end