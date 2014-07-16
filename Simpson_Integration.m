function simp=Simpson_Integration(x,y)
    h = abs(x(2)-x(1));
    n = length(x)-1;
    simp = h*(sum(2*(1+rem(1:n-1,2)).*y(2:n))+y(1)+y(n+1))/3;
end