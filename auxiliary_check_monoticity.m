function monotone=auxiliary_check_monoticity(x)
    L = length(x);
    monotone = true;
    for i=1:L-1
        if x(i) < x(i+1)
            monotone = false;
        end
    end
    if monotone
        return
    else
        monotone = true;
        for i=1:L-1
            if x(i) > x(i+1)
                monotone = false;
            end
        end
    end
end