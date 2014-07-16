function [x_c,y_c,x_sep]=flatten_phase_curve(wavelength,phase,idxs)

L = length(idxs);
x = zeros(1,L);
y = zeros(1,L);
for i=1:L
    idx = find_closest_idx(wavelength,idxs(i));
    if idx==-1
        error('Index not in range!')
    end
    x(i) = wavelength(idx);
    y(i) = phase(idx);
end
x_sep = x(2:end-1);
[x_c,y_c] = calc_catmull_curve(x,y);
end