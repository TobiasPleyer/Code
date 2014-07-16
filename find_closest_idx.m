function idx=find_closest_idx(array, value)
    temp = abs(array-value);
    [temp,idx] = min(temp);
end