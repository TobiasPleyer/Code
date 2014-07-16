function [x_out,y_out]=add_bumb(x_c,y_c,x_o,y_o,idxs,bumbs,center,type,width,height)
    
    x_out = [];
    y_out = [];
    if length(bumbs) ~= length(idxs)-1
        error('The amaount of original bumbs does not correspond to provided list.\nExpected %g',length(idxs)-1)
    end
    original_idxs = arrayfun(@(x) find_closest_idx(x_o,x),idxs);
    catmull_idxs  = arrayfun(@(x) find_closest_idx(x_c,x),idxs);
    for i=1:length(bumbs)
        if bumbs(i)==1
            x_out = [x_out x_o(original_idxs(i):original_idxs(i+1))];
            y_out = [y_out y_o(original_idxs(i):original_idxs(i+1))];
        else
            x_out = [x_out x_c(catmull_idxs(i):catmull_idxs(i+1))];
            y_out = [y_out y_c(catmull_idxs(i):catmull_idxs(i+1))];
        end
    end
    for i=1:length(center)
        half = width(i)/2;
        fprintf('half: %g\n',half)
        idx_lower  = find_closest_idx(x_out,center(i)-half);
        fprintf('idx_lower: %g\n',idx_lower)
        idx_upper  = find_closest_idx(x_out,center(i)+half);
        fprintf('idx_upper: %g\n',idx_upper)
        idx_middle = find_closest_idx(x_out,center(i));
        fprintf('idx_middle: %g\n',idx_middle)
        x_low = x_out(idx_lower);
        fprintf('x_low: %g\n',x_low)
        x_upp = x_out(idx_upper);
        fprintf('x_upp: %g\n',x_upp)
        x_mid = x_out(idx_middle);
        fprintf('x_mid: %g\n',x_mid)
        y_low = y_out(idx_lower);
        fprintf('y_low: %g\n',y_low)
        y_upp = y_out(idx_upper) - y_low;
        fprintf('y_upp: %g\n',y_upp)
        y_mid = height(i);
        fprintf('y_mid: %g\n',y_mid)
        x_mid = x_mid - x_low;
        x_upp = x_upp - x_low;
        coeffs = polyfit([0 x_mid x_upp],[0 y_mid y_upp],2);
        fprintf('a: %g, b: %g, c: %g\n',coeffs(1),coeffs(2),coeffs(3))
        u = x_out(idx_lower:idx_upper) - x_out(idx_lower);
        v = polyval(coeffs,u);
        y_out(idx_lower:idx_upper) = v + y_low;
    end
end