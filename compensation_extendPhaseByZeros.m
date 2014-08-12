function [w_out,I_out,p_out,l_out]=compensation_extendPhaseByZeros(w,I,p,N)
    p_out      = [p zeros(1,N)];
    I_out      = [I zeros(1,N)];
    w_spacing  = w(2)-w(1);
    last       = (w(end)-w(1))/w_spacing;
    w_extended = (last+1):1:(last+N);
    w_extended = w_extended.*w_spacing;
    w_extended = w_extended + w(1);
    w_out      = [w w_extended];
    
    c          = 299792458;
    l_out      = 2*pi*c ./ (w_out*1e15);
end