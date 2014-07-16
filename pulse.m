%% return pulse
function  E = pulse(t, T, lambda0)
    % t in fs
    % T in fs
    % lambda0 in nm
    
    c = 299792458;
    omega0 = 2*pi*c/(lambda0*1e-9)*1e-15;

    E = exp(-2*log(2)*(t./T).^2 - 1i.*(omega0.*t));   
end
