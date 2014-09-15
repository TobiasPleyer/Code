function isIn=auxiliary_isInArray(arr,item)
    isIn = false;
    for i=1:length(arr)
        if arr(i)==item
            isIn = true;
            break;
        end
    end
end