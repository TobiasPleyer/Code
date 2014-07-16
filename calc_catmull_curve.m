function [x_out,y_out]=calc_catmull_curve(x_in,y_in)
    L = length(x_in);
    u = 0:0.01:1;
    l = length(u);
    x_out = zeros(1,(L-3)*l);
    y_out = zeros(1,(L-3)*l);
    %Catmull-Rom curves
    car1 = @(x) -x.^3./2 + x.^2 - x./2;
    car2 = @(x) 3*x.^3./2 - 5*x.^2./2 + 1;
    car3 = @(x) -3*x.^3./2 + 2.*x.^2 + x./2;
    car4 = @(x) x.^3./2 - x.^2./2;
    
    for i=2:L-2
        x_c = x_in(i-1)*car1(u) + x_in(i)*car2(u) + x_in(i+1)*car3(u) + x_in(i+2)*car4(u);
        y_c = y_in(i-1)*car1(u) + y_in(i)*car2(u) + y_in(i+1)*car3(u) + y_in(i+2)*car4(u);
        x_out((i-2)*l+1:(i-1)*l) = x_c;
        y_out((i-2)*l+1:(i-1)*l) = y_c;
    end
end