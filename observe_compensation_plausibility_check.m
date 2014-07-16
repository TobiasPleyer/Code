set(0,'DefaultFigureWindowStyle','docked')

fourier_file = '../Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Et-fourier.dat';
Et_fourier   = dlmread(fourier_file);
t_fourier = Et_fourier(:,1);
I_fourier = Et_fourier(:,2);
int_fourier  = -1*trapz(t_fourier,I_fourier);
clear Et_fourier

filename               = '../Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck-original.dat';
Sk                     = dlmread(filename);
original_wavelength    = Sk(:,1)' .* 1e-9;
original_intensity     = Sk(:,2)';
original_phase         = Sk(:,3)';
clear Sk
%%


%% Physical units, transformations and physical calculations

c                          = 299792458;
original_omega             = 2*pi*c ./ original_wavelength;
indx_lower                 = find_closest_idx(original_omega,2*pi*c/1050e-9);
indx_higher                = find_closest_idx(original_omega,2*pi*c/1018e-9);
cropped_original_omega     = original_omega(indx_lower:indx_higher);
cropped_original_phase     = original_phase(indx_lower:indx_higher);
cropped_original_intensity = original_intensity(indx_lower:indx_higher);
%%


%% Create the test omegas

fun                            = @(x) 250*(x*1e-15-1.79)-8;
plot2_phase                    = arrayfun(fun,cropped_original_omega);
phase2                         = zeros(size(original_phase));
phase2(1:indx_lower-1)         = original_phase(1:indx_lower-1);
phase2(indx_lower:indx_higher) = plot2_phase(:);
phase2(indx_higher+1:end)      = original_phase(indx_higher+1:end);

fun                            = @(x) -25000*(x*1e-15-1.8225)^2+8;
plot3_phase                    = arrayfun(fun,cropped_original_omega);
phase3                         = zeros(size(original_phase));
phase3(1:indx_lower-1)         = original_phase(1:indx_lower-1);
phase3(indx_lower:indx_higher) = plot3_phase(:);
phase3(indx_higher+1:end)      = original_phase(indx_higher+1:end);

fun                            = @(x) 0.8*(5*sin(100*(x*1e-15-1.8225))+5*cos(300*(x*1e-15-1.8225))+5*(sin(900*(x*1e-15-1.8225)))^2-(500*(x*1e-15-1.79)-8)+5*sin(10000*(x*1e-15-1.8225))+4);
plot4_phase                    = arrayfun(fun,cropped_original_omega);
phase4                         = zeros(size(original_phase));
phase4(1:indx_lower-1)         = original_phase(1:indx_lower-1);
phase4(indx_lower:indx_higher) = plot4_phase(:);
phase4(indx_higher+1:end)      = original_phase(indx_higher+1:end);
%%


%% Fourier back transformation


Sk_cplx = sqrt(original_intensity) .* exp(1i*original_phase);
[t,E] = Speck_Fourier(original_wavelength,Sk_cplx);
t1 = t*1e15;
I1 = abs(E).^2;
I1 = I1./max(I1);
int = trapz(t1,I1);
integral_scaling_factor = int_fourier / int;
I1 = I1 .* integral_scaling_factor;
                            
Sk_cplx = sqrt(original_intensity) .* exp(1i*phase2);
[t,E] = Speck_Fourier(original_wavelength,Sk_cplx);
t2 = t*1e15;
I2 = abs(E).^2;
I2 = I2./max(I2);
int = trapz(t2,I2);
integral_scaling_factor = int_fourier / int;
I2 = I2 .* integral_scaling_factor;

Sk_cplx = sqrt(original_intensity) .* exp(1i*phase3);
[t,E] = Speck_Fourier(original_wavelength,Sk_cplx);
t3 = t*1e15;
I3 = abs(E).^2;
I3 = I3./max(I3);
int = trapz(t3,I3);
integral_scaling_factor = int_fourier / int;
I3 = I3 .* integral_scaling_factor;

Sk_cplx = sqrt(original_intensity) .* exp(1i*phase4);
[t,E] = Speck_Fourier(original_wavelength,Sk_cplx);
t4 = t*1e15;
I4 = abs(E).^2;
I4 = I4./max(I4);
int = trapz(t4,I4);
integral_scaling_factor = int_fourier / int;
I4 = I4 .* integral_scaling_factor;
%%


%% Figure part
%-------------------------Figure(1)----------------------------------------
%--------------------------------------------------------------------------
figure(1)
%---------------------------------subplot(4,2,1)---------------------------
subplot(4,2,1)
plot(cropped_original_omega*1e-15,(cropped_original_intensity*25)-15,'color','green')
hold on
plot(cropped_original_omega*1e-15,cropped_original_phase,'color','blue')
plot(cropped_original_omega*1e-15,zeros(size(cropped_original_phase)),'color','black')
hold off
axis([1.79 1.855 -15 15])
title('Original Phase')
xlabel(sprintf('angular frequency in 10^{15}s^{-1}'))
legend('spectrum','phase','zero line','Location','NorthWest')
%---------------------------------subplot(4,2,2)---------------------------
subplot(4,2,2)
plot(t_fourier,I_fourier,'color','blue')
hold on
plot(t1,I1,'color','red')
hold off
axis([-2000 2000 0 1])
title('Time Domain')
xlabel('Time in fs')
legend('fourier','result')
%---------------------------------subplot(4,2,3)---------------------------
subplot(4,2,3)
plot(cropped_original_omega*1e-15,(cropped_original_intensity*25)-15,'color','green')
hold on
plot(cropped_original_omega*1e-15,plot2_phase,'color','blue')
plot(cropped_original_omega*1e-15,zeros(size(cropped_original_phase)),'color','black')
hold off
axis([1.79 1.855 -15 15])
title('Linear Phase')
xlabel(sprintf('angular frequency in 10^{15}s^{-1}'))
legend('spectrum','phase','zero line','Location','NorthWest')
%---------------------------------subplot(4,2,4)---------------------------
subplot(4,2,4)
plot(t_fourier,I_fourier,'color','blue')
hold on
plot(t2,I2,'color','red')
hold off
axis([-2000 2000 0 1])
title('Time Domain')
xlabel('Time in fs')
legend('fourier','result')
%---------------------------------subplot(4,2,5)---------------------------
subplot(4,2,5)
plot(cropped_original_omega*1e-15,(cropped_original_intensity*25)-15,'color','green')
hold on
plot(cropped_original_omega*1e-15,plot3_phase,'color','blue')
plot(cropped_original_omega*1e-15,zeros(size(cropped_original_phase)),'color','black')
hold off
axis([1.79 1.855 -15 15])
title('Quadratic Phase')
xlabel(sprintf('angular frequency in 10^{15}s^{-1}'))
legend('spectrum','phase','zero line')
%---------------------------------subplot(4,2,6)---------------------------
subplot(4,2,6)
plot(t_fourier,I_fourier,'color','blue')
hold on
plot(t3,I3,'color','red')
hold off
axis([-2000 2000 0 1])
title('Time Domain')
xlabel('Time in fs')
legend('fourier','result')
%---------------------------------subplot(4,2,7)---------------------------
subplot(4,2,7)
plot(cropped_original_omega*1e-15,(cropped_original_intensity*25)-15,'color','green')
hold on
plot(cropped_original_omega*1e-15,plot4_phase,'color','blue')
plot(cropped_original_omega*1e-15,zeros(size(cropped_original_phase)),'color','black')
hold off
axis([1.79 1.855 -15 15])
title('Complex Phase')
xlabel(sprintf('angular frequency in 10^{15}s^{-1}'))
legend('spectrum','phase','zero line')
%---------------------------------subplot(4,2,8)---------------------------
subplot(4,2,8)
plot(t_fourier,I_fourier,'color','blue')
hold on
plot(t4,I4,'color','red')
hold off
axis([-2000 2000 0 1])
title('Time Domain')
xlabel('Time in fs')
legend('fourier','result')