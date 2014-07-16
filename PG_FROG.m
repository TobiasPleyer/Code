function [tau_array, lambda_array, spectro, type] = PG_FROG(E_vec, tau_range, tbuf,...
                                        lambda0, l_min, l_max, gridSize)
    %calculated the spectrogram of a PG-FROG
    
    type = 'PG-FROG';
    
    c = 299792458e-6; %nm/fs
    
    f_min = c.*(1/(1.01.*l_max));
    f_max = c.*(1/(0.99.*l_min));
    
    N_tau = gridSize;
    N_f =  gridSize;
    
    tau_array = linspace(-tau_range, tau_range, N_tau);
    f_array = linspace(f_min, f_max, N_f);
    lambda_array = c./f_array;  %nm
    
    spectrogram = zeros(N_f,N_tau);
    S_lambda = zeros(N_f,N_tau);
    
    t_max = tbuf+tau_range;
    dt = 1/(2.*f_max);
    N_t = 2^nextpow2(2.*t_max/dt);
    dt = 2*t_max/(N_t-1);
    f_fft = 1/(2*dt).*linspace(0,1,N_t/2+1);
    
    t = linspace(-t_max,t_max,N_t);
    E = interp1(E_vec(:,1), E_vec(:,2), t);
    Phi = interp1(E_vec(:,1), E_vec(:,3), t);
    E = sqrt(E).*exp(1i.*Phi).*exp(-1i.*2*pi.*c/lambda0.*t);
    E(isnan(E)) = 0;
    
    %[~, E] = GDDpulse(100, 1030, 5e4, t_max, N_t);
    
    for n=1:N_tau
        tau = tau_array(n);
        E2 = interp1(t,E,t-tau,'nearest');
        E2(isnan(E2)) = 0;
        %E2 = pulse(t-tau, 200, 1030);
        Y = E.*abs(E2).^2;
        
        S = abs(ifft(fftshift(Y)));
        
        S = S(1:N_t/2+1).^2;
        spectrogram(:,n) = interp1(f_fft, S, f_array);
        S_lambda(:,n) = c./lambda_array'.^2.*spectrogram(:,n);
    end
    spectro = S_lambda;
end