%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Written by Tobias Pleyer 2014
%       Version 1.0
%
%       This script tries to evaluate if the measured FROG pahse is given
%       with the correct sign or if it has to be flipped.
%
%       VERSION
%       1.0 Multiply the effect of graating to phase and observe if a
%           compression is achieved. (2014/06/05)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

set(0,'DefaultFigureWindowStyle','docked')

%% This is the part where all filenames and directories are defined

parent    = '../Daten/';
filebase  = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';

%%


%% Open the files and get the original data
    
    % Wavelength based field
    filename_Speck          = sprintf('%s%s.Speck.dat',parent,filebase);
    Sk                      = dlmread(filename_Speck);
    original_wavelength     = Sk(:,1)' .* 1e-9;
    original_intensity      = Sk(:,2)';
    original_phase          = Sk(:,3)';
    
    % Calculate the Fourier limit of the pulse
    Sk_cplx                 = sqrt(original_intensity);
    [t_fourier,E_fourier]   = Speck_Fourier(original_wavelength,Sk_cplx);
    t_fourier               = t_fourier*1e15;
    I_fourier               = abs(E_fourier).^2;
    I_fourier               = I_fourier./max(I_fourier);
    Int_fourier             = trapz(t_fourier,I_fourier);
    
    % Time based field, phase compensated
    filename_Ek             = sprintf('%s%s.Ek.dat',parent,filebase);
    Ek                      = dlmread(filename_Ek);
    t_Ek                    = Ek(:,1);
    I_Ek                    = Ek(:,2);
    Int_Ek                  = abs(trapz(t_Ek,I_Ek));
    integral_scaling_factor = Int_fourier / Int_Ek;   
%%


%% Physical units, transformations and physical calculations

    c                          = 299792458;
    omega0                     = 2*pi*c ./ 1.030e-6;
    original_omega             = 2*pi*c ./ (original_wavelength); % needed to map the correct phase to the corresponding frequency, hence the fliplr
    indx_lower                 = find_closest_idx(original_omega,1.80e15);
    indx_higher                = find_closest_idx(original_omega,1.855e15);
    cropped_original_omega     = original_omega(indx_lower:indx_higher);
    cropped_original_phase     = original_phase(indx_lower:indx_higher);
    cropped_original_intensity = original_intensity(indx_lower:indx_higher);
%%


%% Calculate the contribution of the grating

    N = 300/0.001; % 3e5
    beta0 = asin(N*1.030e-6 - sin(2*pi/360*30)); % estimated angle of incidence 30°
    D = 0.12; % 12cm
    D2 = -2.*D./c.*1./original_omega.*((2*pi*c*N)./(original_omega.*cos(beta0))).^2;
    phi = 0.5.*D2.*(original_omega-omega0).^2;
%%

    
%% Fourier back integration part

    Sk_cplx = sqrt(original_intensity) .* exp(1i*original_phase) .* exp(1i*phi);
    [t_positive,E] = Speck_Fourier(original_wavelength,Sk_cplx);
    t_positive = t_positive*1e15;
    I_positive = abs(E).^2;
    I_positive = I_positive./max(I_positive);
    Int_Ek = trapz(t_positive,I_positive);
    integral_scaling_factor = Int_fourier / Int_Ek;
    %I = I .* integral_scaling_factor;
    Sk_cplx = sqrt(original_intensity) .* exp(-1i*original_phase) .* exp(1i*phi);
    [t_negative,E] = Speck_Fourier(original_wavelength,Sk_cplx);
    t_negative = t_negative*1e15;
    I_negative = abs(E).^2;
    I_negative = I_negative./max(I_negative);
    Int_Ek = trapz(t_negative,I_negative);
    integral_scaling_factor = Int_fourier / Int_Ek;
%%


%% Plotting part

figure(1)
plot(t_fourier,I_fourier,'color','blue')
hold on
plot(t_Ek,I_Ek,'color','green')
plot(t_positive,I_positive,'color','red')
plot(t_negative,I_negative,'color','black')
hold off
legend('fourier limit','original','normal compression','flipped compression')
axis([-2500 4000 0 1])
text(2000,0.8,sprintf('Calculate the contribution of the grating\nN = 300/0.001;\nbeta0 = asin(N*1.030e-6 - sin(2*pi/360*30));\nD = 0.12;\nD2 = -2.*D./c.*1./original_omega.*((2*pi*c*N)./(original_omega.*cos(beta0))).^2;\nphi = 0.5.*D2.*(original_omega-omega0).^2;\n'),'FontSize',8)