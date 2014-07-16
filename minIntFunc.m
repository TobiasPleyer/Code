function F=minIntFunc(x,omega,phase,gdd)
    % To avoid dimension problems
    omega = omega(:);
    phase = phase(:);
    gdd   = gdd(:);
    % Integration part
    I = cumtrapz(omega,gdd);
    I = cumtrapz(omega,I+x(1));
    I = I*1e30+x(2);
    % Optimization part
    F = sum(abs(phase-I));
end