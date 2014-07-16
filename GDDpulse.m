
%% return pulse
function [t1,E] = GDDpulse(T, lambda0, GDD, t_max, N_t)
    % T in fs
    % lambda0 in nm
    % GDD in fs^2
    % N number of points in time domain
    % WIN width of time window in multiples of T
    
    c = 299792458;
    f_c = c/(lambda0*1e-9)*1e-15;
    omega0 = 2*pi*f_c;

    dt = 2*t_max/(N_t-1);
    t1 = linspace(-t_max,t_max,N_t);
    E1 = exp(-2*log(2)*(t1./T).^2 - 1i.*(omega0.*t1));   
    
    %f_c = omega0/2/pi;
    S = ifft(ifftshift(E1));
    f = 1/(2*dt).*linspace(0,1,N_t/2+1);
    
    H_ef = exp(-1i.*2.*pi^2.*GDD.*(f-f_c).^2);
    S(1:N_t/2+1) = S(1:N_t/2+1).*H_ef;
    
    E = fftshift(fft(S));
   % plotyy(f,abs(S(1:size(f,2))).^2, f, H_ef)
   % plot(t1,abs(E).^2)
end

