function Y=adjust_height(x,y,lmin,lmax)
    Y = y(:);
    lidx = find_closest_idx(x,lmin); %lower index
    uidx = find_closest_idx(x,lmax); %upper index
    diff = Y(uidx)-Y(lidx);
    Y(uidx:end) = Y(uidx:end) - diff;
    idx_range = lidx:1:uidx;
    xrange = x(idx_range);
    xmin = x(lidx);
    xmax = x(uidx);
    ymin = Y(lidx);
    ymax = Y(uidx);
    yrange = interp1([xmin xmax],[ymin ymax],xrange);
    Y(idx_range) = yrange;
end