function GDD=calc_GDD_from_phase(phase,h)
    GDD = (phase(3:end)+phase(1:end-2)-2*phase(2:end-1))./h^2;
end