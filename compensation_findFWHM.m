function [t_min,t_max]=compensation_findFWHM(time,intensity)
    normed_intensity = intensity./max(intensity);
    above = time(normed_intensity>=0.5);
    t_min = min(above);
    t_max = max(above);
end