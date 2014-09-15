function [t_min,t_max,I_minmax]=auxiliary_findFWHM(time,intensity)
    M = max(intensity);
    normed_intensity = intensity./M;
    above = time(normed_intensity>=0.5);
    t_min = min(above);
    t_max = max(above);
    I_minmax = 0.5 * M;
end