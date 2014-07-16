%% ########################################################################
%
%  
%  Version: 1.0
%
%  INPUT:
%
%  OUTPUT:
%
%  HISTORY:
%      v1.0: First runable state (16.05.2014)
%
%% ########################################################################

%% INITILIZATION

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
%format compact;
set(0,'DefaultFigureWindowStyle','docked')

% Change the current folder to the folder of this m-file.
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end
%%
%% This is the part where all filenames and directories are defined

parent    = '../Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/';
filebase  = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
filebase2 = '26_AIR_FROG_31-W_RTT=7-us_Ip=12-A-bin';
endings   = {'center-flattened', ...
             'only-third-bumb', ...
             'smooth', ...
             'wings-flattened', ...
             'catmull', ...
             'smooth-with-randn-pm05', ...
             'smooth-with-randn-pm005', ...
             'original-with-randn-pm05', ...
             'original-with-randn-pm005', ...
             'adjust-height-10241-10248', ...
             'adjust-height-10268-10275'};             

%%


%% Fourier limit stuff which is independent of the loop

    fourier_file = sprintf('%s%s.Et-%s.dat',parent,filebase,'fourier');
    Et_fourier   = dlmread(fourier_file);
    int_fourier  = abs(trapz(Et_fourier(:,1),Et_fourier(:,2)));
    trebino_file = sprintf('%s%s.Ek-%s.dat',parent,filebase,'trebino');
    Ek_trebino   = dlmread(trebino_file);
    
    t_fourier = Et_fourier(:,1);
    I_fourier = Et_fourier(:,2);
%%

for i=3%:length(endings)
    fprintf('Start...%s\n',endings{i})
%%


%% Open the files and get the original data
    
    % Time based field, phase compensated
    filename_Et             = sprintf('%s%s.Et-%s.dat',parent,filebase,endings{i});
    Et                      = dlmread(filename_Et);
    t                       = Et(:,1);
    I_smoothed              = Et(:,2);
    int                     = abs(trapz(t,I_smoothed));
    integral_scaling_factor = int_fourier / int;
    I_smoothed              = I_smoothed .* integral_scaling_factor;
    
    % Time based field, phase uncompensated
    filename_Ek = sprintf('%s%s.Ek-%s.dat',parent,filebase,endings{i});
    Ek = dlmread(filename_Ek);
    
    % Wavelength based field
    name                   = sprintf('../Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck-%s.dat',endings{i});
    Sk                     = dlmread(name);
    original_wavelength    = fliplr(Sk(:,1)' .* 1e-9);
    original_phase         = fliplr(Sk(:,3)' + Sk(:,4)'); % Because SK(:,3)=new_phase and Sk(:,4)=old_phase-new_phase
    interpolated_new_phase = fliplr(Sk(:,3)');
    original_intensity     = fliplr(Sk(:,2)');
%%


%% Physical units, transformations and physical calculations

    c                          = 299792458;
    original_omega             = 2*pi*c ./ (original_wavelength); % needed to map the correct phase to the corresponding frequency, hence the fliplr
    indx_lower                 = find_closest_idx(original_omega,1.80e15);
    indx_higher                = find_closest_idx(original_omega,1.855e15);
    cropped_original_omega     = original_omega(indx_lower:indx_higher);
    cropped_original_phase     = original_phase(indx_lower:indx_higher);
    cropped_original_intensity = original_intensity(indx_lower:indx_higher);
    cropped_new_phase          = interpolated_new_phase(indx_lower:indx_higher); 
%%


%% Differentiation part

    increment              = 0.5e11;
    start                  = cropped_original_omega(1);
    ende                   = cropped_original_omega(end);
    vq                     = start:increment:ende;
    interpolated_phase     = interp1(cropped_original_omega,cropped_original_phase,vq,'spline');
    interpolated_new_phase = interp1(cropped_original_omega,cropped_new_phase,vq,'spline');
    GDD_original_phase     = calc_GDD_from_phase(interpolated_phase,increment);
    GDD_new_phase          = calc_GDD_from_phase(interpolated_new_phase,increment);
    smoothed               = smooth(GDD_new_phase,80,'sgolay',4);
    super_smoothed         = smooth(smoothed,50);
    interpolated_new_phase = interpolated_new_phase(2:end-1);
    interpolated_omega     = vq(2:end-1) .* 1e-15;
%%


%% Low pass filtering
    %-------------------------make preparations for a lowpass filter-------
    sampling_freq = 1/(interpolated_omega(2)-interpolated_omega(1)); % 0.2THz here
    order         = 2;
    cut_off_freq  = 50;
    peak_to_peak_dB = 0.5;      % This is the Matlab doc recommend first guess
    stopband_atten  = 20;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A] = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
    N = 5000;                   % number of points used to create the filter
                                % show the figures of merit for the filter
                                % figure()
                                % freqz(B,A,N,sampling_freq)
    %-------------------------filter the data------------------------------
                                % filtfilt has the big advantage over filter that it operates at a zero phase
                                % shift, thus allowing to obtain a direct correspondance to the filtered data
    filtered_super_smoothed = filtfilt(B,A,super_smoothed);
    %-------------------------same procedure for the original--------------
    sampling_freq = 1/(interpolated_omega(2)-interpolated_omega(1)); % 0.2THz here
    order         = 2;
    cut_off_freq  = 60;
    peak_to_peak_dB = 0.8;      % This is the Matlab doc recommend first guess
    stopband_atten  = 50;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A] = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
    filtered_original = filtfilt(B,A,GDD_original_phase);
%%


%% Integration part

    sol = fminsearch(@(x)minIntFunc(x,interpolated_omega,interpolated_phase(2:end-1),smoothed),[1e-28,3.0]);
    smoothed_GDD_integrated = cumtrapz(interpolated_omega,smoothed);
    smoothed_GDD_integrated = cumtrapz(interpolated_omega,smoothed_GDD_integrated+sol(1));
    smoothed_GDD_integrated = smoothed_GDD_integrated*1e30+sol(2);
    
    sol = fminsearch(@(x)minIntFunc(x,interpolated_omega,interpolated_phase(2:end-1),super_smoothed),[1e-28,3.0]);
    super_smoothed_integrated = cumtrapz(interpolated_omega,super_smoothed);
    super_smoothed_integrated = cumtrapz(interpolated_omega,super_smoothed_integrated+sol(1));
    super_smoothed_integrated = super_smoothed_integrated*1e30+sol(2);
    
    sol = fminsearch(@(x)minIntFunc(x,interpolated_omega,interpolated_phase(2:end-1),filtered_super_smoothed),[1e-29,3.0]);
    filtered_super_smoothed_integrated = cumtrapz(interpolated_omega,filtered_super_smoothed);
    filtered_super_smoothed_integrated = cumtrapz(interpolated_omega,filtered_super_smoothed_integrated+sol(1));
    filtered_super_smoothed_integrated = filtered_super_smoothed_integrated*1e30+sol(2);
    
    sol = fminsearch(@(x)minIntFunc(x,interpolated_omega,interpolated_phase(2:end-1),filtered_original),[1e-28,3.0]);
    filtered_original_integrated = cumtrapz(interpolated_omega,filtered_original);
    filtered_original_integrated = cumtrapz(interpolated_omega,filtered_original_integrated+sol(1));
    filtered_original_integrated = filtered_original_integrated*1e30+sol(2);
%%


%% Resampling
    
                                % Interpolate it and bring it on the same sampling rate as the original
                                % phase
                                % awkward; the reason for this indexing lies in the fact that diff reduces the length of an
                                % array by 1, as in the calculation of the GDD. See also the indexing of omega
    truncated_super_smoothed_integrated           = interp1(interpolated_omega,super_smoothed_integrated,cropped_original_omega(2:end-1)*1e-15);
                                
    truncated_filtered_super_smoothed_integrated  = interp1(interpolated_omega,filtered_super_smoothed_integrated,cropped_original_omega(2:end-1)*1e-15);
    
    truncated_filtered_original_integrated        = interp1(interpolated_omega,filtered_original_integrated,cropped_original_omega(2:end-1)*1e-15);
%%


%% Combine the new phase curve part with original data

                                % This indexing is more than ugly but for now the best option to reverse the phase shift induced by the filtering
    truncated_filtered_super_smoothed_integrated_patched                              = zeros(size(original_phase));
    truncated_filtered_super_smoothed_integrated_patched(1:indx_lower)                = original_phase(1:indx_lower);
    truncated_filtered_super_smoothed_integrated_patched(indx_lower+1:indx_higher-1)  = truncated_filtered_super_smoothed_integrated;
    truncated_filtered_super_smoothed_integrated_patched(indx_higher:end)             = original_phase(indx_higher:end);
    cropped_truncated_filtered_super_smoothed_integrated_patched                      = truncated_filtered_super_smoothed_integrated_patched(indx_lower:indx_higher);
    
    truncated_super_smoothed_integrated_patched                             = zeros(size(original_phase));
    truncated_super_smoothed_integrated_patched(1:indx_lower)               = original_phase(1:indx_lower);
    truncated_super_smoothed_integrated_patched(indx_lower+1:indx_higher-1) = truncated_filtered_super_smoothed_integrated;
    truncated_super_smoothed_integrated_patched(indx_higher:end)            = original_phase(indx_higher:end);
    cropped_truncated_super_smoothed_integrated_patched                     = truncated_filtered_super_smoothed_integrated_patched(indx_lower:indx_higher);
    
    truncated_filtered_original_integrated_patched                             = zeros(size(original_phase));
    truncated_filtered_original_integrated_patched(1:indx_lower)               = original_phase(1:indx_lower);
    truncated_filtered_original_integrated_patched(indx_lower+1:indx_higher-1) = truncated_filtered_original_integrated;
    truncated_filtered_original_integrated_patched(indx_higher:end)            = original_phase(indx_higher:end);
    cropped_truncated_filtered_original_integrated_patched                     = truncated_filtered_original_integrated_patched(indx_lower:indx_higher);
%% Fourier back transformation

                                % Compare the effect of compensating the origninal phase by the
                                % modified phase, brought to the same sampling rate
    diff_original_to_super_smooth = original_phase - truncated_super_smoothed_integrated_patched;
    cropped_diff_original_to_super_smooth = diff_original_to_super_smooth(indx_lower:indx_higher);
    Sk_cplx = sqrt(original_intensity) .* exp(1i*diff_original_to_super_smooth);
    [t,E] = Speck_Fourier(original_wavelength,Sk_cplx);
    t = t*1e15;
    I = abs(E).^2;
    I = I./max(I);
    int = trapz(t,I);
    integral_scaling_factor = int_fourier / int;
    I = I .* integral_scaling_factor;
                                % Same as above, only that the compensating phase is the phase we
                                % obtained by filtering the GDD
    diff_original_to_filtered_original = original_phase - truncated_filtered_original_integrated_patched;
    cropped_diff_original_to_filtered_original = diff_original_to_filtered_original(indx_lower:indx_higher);
    Sk_cplx_filtered = sqrt(original_intensity) .* exp(1i*diff_original_to_filtered_original);
    [t_filtered,E_filtered] = Speck_Fourier(original_wavelength,Sk_cplx_filtered);
    t_filtered = t_filtered*1e15;
    I_filtered = abs(E_filtered).^2;
    I_filtered = I_filtered./max(I_filtered);
    int_filtered = trapz(t_filtered,I_filtered);
    factor_filtered = int_fourier / int_filtered;
    I_filtered = I_filtered .* factor_filtered;
%%


%% Figure part
    %-------------------------Figure(1)------------------------------------
    %----------------------------------------------------------------------
    figure(1)
    %---------------------------------subplot(2,2,1)-----------------------
    subplot(2,2,1)
    plot(t_fourier,I_fourier,'color','red');
    hold 'on'
    plot(t,I_smoothed,'color','green');
    % We are ready to add the additional plots to the Fourier limit
    plot(t,I,'color','magenta')
    plot(t_filtered,I_filtered,'color','black')
    title('Comparison between the compensated pulses and the fourier limit')
    text(500,0.75,sprintf('Curves have the\nsame integral'),'FontSize',5,'BackgroundColor',[.7 .9 .7])
    axis([-2000 2000 0 1])
    legend('fourier limit','smooth phase','super smooth','filtered original')
    hold 'off'
    %---------------------------------subplot(2,2,2)-----------------------
    subplot(2,2,2)
    plot(interpolated_omega,GDD_original_phase.*1e30, 'r')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    title('Comparison of all GDD curves')
    axis([1.80 1.855 -1.5e6 1.5e6])
    hold 'on'
    plot(interpolated_omega,GDD_new_phase.*1e30)
    plot(interpolated_omega,smoothed.*1e30,'color','green')
    plot(interpolated_omega,super_smoothed*1e30,'color','magenta','LineWidth',3)
    plot(interpolated_omega,filtered_original*1e30,'color','black','LineWidth',3)
    hold 'off'
    legend('original','from smoothed phase','smoothed GDD from smoothed phase','super smooth','filtered original','Location','NorthOutside')
    %---------------------------------subplot(2,2,3)-----------------------
    subplot(2,2,3)
     plot(interpolated_omega,super_smoothed*1e30,'color','green','LineWidth',3)
    hold on
    plot(interpolated_omega,filtered_super_smoothed*1e30,'r')
    plot(interpolated_omega,filtered_original*1e30,'color','black')
    hold off
    legend('super smoothed','super smoothed filtered','original filtered','Location','NorthWest')
    title('Comparison of the best options for a simple GDD')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    xlim([1.80 1.855])
    %---------------------------------subplot(2,2,4)-----------------------
    subplot(2,2,4)
    plot(cropped_original_omega*1e-15,cropped_original_phase,'color','cyan')
    hold on
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_super_smoothed_integrated_patched,'color','magenta','LineWidth',3)
    plot(cropped_original_omega*1e-15,cropped_truncated_super_smoothed_integrated_patched,'color','black')
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_original_integrated_patched,'color','red')
    plot(cropped_original_omega*1e-15,cropped_diff_original_to_super_smooth,'color','black')
    plot(cropped_original_omega*1e-15,cropped_diff_original_to_filtered_original,'color','red') 
    plot(cropped_original_omega.*1e-15,(cropped_original_intensity .* 22)-13,'color','green')
    % Show the zero line
    plot(cropped_original_omega*1e-15,zeros(size(cropped_diff_original_to_filtered_original)),'color','black')
    hold off
    axis([1.80 1.855 -13 11])
    legend('original','filtered super smoothed','super smoothed','filtered original','difference to super smoothed','difference to filtered original','spectrum','zero line','Location','EastOutside')
    title('Direct comparison between the original phase curve and the integrated alternatives (integration constants adapted)')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel('Phase')
    %-------------------------Figure(2)------------------------------------
    %----------------------------------------------------------------------
    figure(2)
    plot(interpolated_omega,GDD_original_phase.*1e30,'color','red')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    title('Comparison of all three GDD curves')
    %axis([1.80 1.855 -1.5e6 1.5e6])
    hold 'on'
    plot(interpolated_omega,GDD_new_phase.*1e30)
    plot(interpolated_omega,smoothed.*1e30,'color','green')
    plot(interpolated_omega,super_smoothed*1e30,'color','magenta','LineWidth',3)
    plot(interpolated_omega,filtered_original*1e30,'color','black','LineWidth',3)
    hold 'off'
    legend('original','from smoothed phase','smoothed GDD from smoothed phase','super smooth','filtered original')
    %-------------------------Figure(3)------------------------------------
    %----------------------------------------------------------------------
    figure(3)
    plot(interpolated_omega,super_smoothed*1e30,'color','green','LineWidth',3)
    hold on
    plot(interpolated_omega,filtered_super_smoothed*1e30,'r')
    plot(interpolated_omega,filtered_original*1e30,'color','black')
    hold off
    legend('super smoothed','super smoothed filtered','original filtered')
    title('Comparison of the best options for a simple GDD')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    xlim([1.80 1.855])
    %-------------------------Figure(4)------------------------------------
    %----------------------------------------------------------------------
    figure(4)
    plot(cropped_original_omega*1e-15,cropped_original_phase,'color','cyan')
    hold on
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_super_smoothed_integrated_patched,'color','magenta','LineWidth',3)
    plot(cropped_original_omega*1e-15,cropped_truncated_super_smoothed_integrated_patched,'color','black')
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_original_integrated_patched,'color','red')
    plot(cropped_original_omega*1e-15,cropped_diff_original_to_super_smooth,'color','black')
    plot(cropped_original_omega*1e-15,cropped_diff_original_to_filtered_original,'color','red') 
    plot(cropped_original_omega.*1e-15,(cropped_original_intensity .* 22)-13,'color','green')
    % Show the zero line
    plot(cropped_original_omega*1e-15,zeros(size(cropped_diff_original_to_filtered_original)),'color','black')
    hold off
    axis([1.80 1.855 -13 11])
    legend('original','filtered super smoothed','super smoothed','filtered original','difference to super smoothed','difference to filtered original','spectrum','zero line')
    title('Direct comparison between the original phase curve and the integrated smoothed GDD (integration constants adapted)')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel('Phase')
%%


%% Clear all unwanted variables

    clear Sk
    clear vq
    %clear increment start ende
%%
    %saveas(gcf,sprintf('Daten/Pictures/sum_up/comparison_to_fourier_%s',endings{i}),'jpg')
end