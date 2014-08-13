function p=compensation_calcPoly3(C,w0)
    a = C(1);
    b = -3*a*w0 + C(2);
    c = 3*a*w0^2 - 2*C(2)*w0;
    d = -a*w0^3 + C(2)*w0^2;
    p = [a b c d];
end