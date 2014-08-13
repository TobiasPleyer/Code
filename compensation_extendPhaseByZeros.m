function [w_ext,I_ext,p_ext,l_ext]=compensation_extendPhaseByZeros(w,I,p,N)
    p_ext      = [p zeros(1,N)];
    I_ext      = [I zeros(1,N)];
    w_spacing  = w(2)-w(1);
    last       = (w(end)-w(1))/w_spacing;
    w_extended = (last+1):1:(last+N);
    w_extended = w_extended.*w_spacing;
    w_extended = w_extended + w(1);
    w_ext      = [w w_extended];
    
    c          = 299792458;
    l_ext      = 1e9.*2*pi*c ./ (w_ext*1e15);
end