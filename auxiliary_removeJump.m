function [out]=auxiliary_removeJump(vector,index,y,n)
    out = vector;
    if n<0
        out(index:end) = out(index:end) + y;
    else
        out(1:index) = out(1:index) + y;
    end
    x1 = out(index);
    x2 = out(index+n);
    inc = (x2-x1)/n;
    if n<0
        a = n:-1;
    else
        a = 1:n;
    end
    for i=a
        out(index+i) = out(index) + i*inc;
    end
end