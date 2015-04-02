function idx=auxiliary_find_closest_idx(array, value)
    temp = abs(array-value);
    [~,idx] = min(temp);
end