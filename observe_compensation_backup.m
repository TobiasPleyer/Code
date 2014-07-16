set(0,'DefaultFigureWindowStyle','docked')
%% This is the part where all filenames and directories are defined

parent    = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/';
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
    int_fourier  = trapz(Et_fourier(:,1),Et_fourier(:,2));
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
    I                       = Et(:,2);
    int                     = trapz(t,I);
    integral_scaling_factor = int_fourier / int;
    I                       = I .* integral_scaling_factor;
    
    % Time based field, phase uncompensated
    filename_Ek = sprintf('%s%s.Ek-%s.dat',parent,filebase,endings{i});
    Ek = dlmread(filename_Ek);
    
    % Wavelength based field
    name                   = sprintf('Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck-%s.dat',endings{i});
    Sk                     = dlmread(name);
    wavelength             = Sk(:,1) .* 1e-9;
    original_phase         = Sk(:,3) + Sk(:,4); % Because SK(:,3)=new_phase and Sk(:,4)=old_phase-new_phase
    interpolated_new_phase = Sk(:,3);
    original_intensity     = Sk(:,2);
    clear Sk
%%


%% Physical units, transformations and physical calculations

    c                          = 299792458;
    original_omega             = fliplr(2*pi*c ./ (wavelength')); % needed to map the correct phase to the corresponding frequency, hence the fliplr
    indx_lower                 = find_closest_idx(original_omega,2*pi*c/1050e-9);
    indx_higher                = find_closest_idx(original_omega,2*pi*c/1018e-9);
    cropped_original_omega     = original_omega(indx_lower:indx_higher);
    cropped_original_phase     = original_phase(indx_lower:indx_higher);
    cropped_original_intensity = original_intensity(indx_lower:indx_higher);
    cropped_new_phase          = interpolated_new_phase(indx_lower:indx_higher); 
%%


%% Differentiation part

    increment              = 0.5e11;
    vq                     = 1.80e15:increment:1.85e15;
    interpolated_phase     = interp1(cropped_original_omega,cropped_original_phase,vq,'spline');
    interpolated_new_phase = interp1(cropped_original_omega,cropped_new_phase,vq,'spline');
    GDD_original_phase     = calc_GDD_from_phase(interpolated_phase,increment);
    GDD_new_phase          = calc_GDD_from_phase(interpolated_new_phase,increment);
    smoothed               = smooth(GDD_new_phase,80,'sgolay',4);
    super_smoothed         = smooth(smoothed,50);
    interpolated_new_phase = interpolated_new_phase(2:end-1);
    interpolated_new_phase = interpolated_new_phase(:);
    omega                  = vq(2:end-1) .* 1e-15;
%%


%% Low pass filtering
    %-------------------------make preparations for a lowpass filter-------
    sampling_freq = 1/(omega(2)-omega(1)); % 0.2THz here
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
    filtered_data = filter(B,A,super_smoothed);
%temp = filtfilt(B,A,super_smoothed);
    %-------------------------same procedure for the original--------------
    sampling_freq = 1/(omega(2)-omega(1)); % 0.2THz here
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
    filtered_original = filter(B,A,GDD_original_phase);
%%


%% Integration part

    a = 2.12e-29;
    b = 4.3;
    smoothed_GDD_integrated = cumtrapz(omega,smoothed);
    smoothed_GDD_integrated = cumtrapz(omega,smoothed_GDD_integrated+a);
    smoothed_GDD_integrated = smoothed_GDD_integrated*1e30+b;
    a = 0.5e-29;
    b = 4.1;
    super_smoothed_integrated = cumtrapz(omega,super_smoothed);
    super_smoothed_integrated = cumtrapz(omega,super_smoothed_integrated+a);
    super_smoothed_integrated = super_smoothed_integrated*1e30+b;
    a = 2.5e-29;
    b = 3.8;
    filtered_data_integrated = cumtrapz(omega,filtered_data);
    filtered_data_integrated = cumtrapz(omega,filtered_data_integrated+a);
    filtered_data_integrated = filtered_data_integrated*1e30+b;
%%


%% Resampling
    
                                % Interpolate it and bring it on the same sampling rate as the original
                                % phase
    truncated_super_smoothed_integrated = zeros(size(original_phase));
    truncated_filtered_data_integrated  = zeros(size(original_phase));
    L                                   = length(super_smoothed_integrated);
    for j=1:length(original_phase)
        idx = find_closest_idx(super_smoothed_integrated,original_phase(j));
        if (idx == 1) || (idx == L)
            truncated_super_smoothed_integrated(j) = original_phase(j);
            truncated_filtered_data_integrated(j)  = original_phase(j);
        else
            truncated_super_smoothed_integrated(j) = super_smoothed_integrated(idx);
            % CAREFUL! MIGHT BE WRONG!
            truncated_filtered_data_integrated(j)  = filtered_data_integrated(idx);
        end
    end
                                % awkward; the reason for this indexing lies in th fact that diff reduces the length of an
                                % array by 1, as in the calculation of the GDD. See also the indexing of omega
    truncated_filtered_data_integrated = interp1(omega,filtered_data_integrated,cropped_original_omega(12:end-1)*1e-15);
%%


%% Combine the new phase curve part with original data

                                % This indexing is more than ugly but for now the best option to reverse the phase shift induced by the filtering
    truncated_filtered_data_integrated_patched                              = zeros(size(original_phase));
    truncated_filtered_data_integrated_patched(1:indx_lower+10)             = original_phase(1:indx_lower+10);
    truncated_filtered_data_integrated_patched(indx_lower+11:indx_higher-7) = truncated_filtered_data_integrated(7:end);
    truncated_filtered_data_integrated_patched(indx_higher-6:end)           = original_phase(indx_higher-6:end);
    cropped_truncated_filtered_data_integrated_patched                      = truncated_filtered_data_integrated_patched(indx_lower:indx_higher);

%% Fourier back transformation

                                % Compare the effect of compensating the origninal phase by the
                                % modified phase, brought to the same sampling rate
    diff_phase = original_phase - truncated_super_smoothed_integrated;
    Sk_cplx = sqrt(original_intensity .* exp(1i*diff_phase));
    [t,E] = Speck_Fourier(wavelength,Sk_cplx);
    t = t*1e15;
    E = abs(E).^2;
    E = E./max(E);
    int = trapz(t,E);
    integral_scaling_factor = int_fourier / int;
    E = E .* integral_scaling_factor;
                                % Same as above, only that the compensating phase is the phase we
                                % obtained by filtering the GDD
    diff_phase_filtered = original_phase(:) - truncated_filtered_data_integrated_patched(:);
    cropped_diff_phase_filtered = diff_phase_filtered(indx_lower:indx_higher);
    Sk_cplx_filtered = sqrt(original_intensity .* exp(1i*diff_phase_filtered));
    [t_filtered,E_filtered] = Speck_Fourier(wavelength,Sk_cplx_filtered);
    t_filtered = t_filtered*1e15;
    E_filtered = abs(E_filtered).^2;
    E_filtered = E_filtered./max(E_filtered);
    int_filtered = trapz(t_filtered,E_filtered);
    factor_filtered = int_fourier / int_filtered;
    E_filtered = E_filtered .* factor_filtered;
%%


%% Figure part
    %-------------------------Figure(1)------------------------------------
    %----------------------------------------------------------------------
    figure(1)
    %---------------------------------subplot(2,2,1)-----------------------
    subplot(2,2,1)
    plot(t_fourier,I_fourier,'color','red');
    hold 'on'
    plot(t,I,'color','green');
    % We are ready to add the additional plots to the Fourier limit
    plot(t,E,'color','magenta')
    plot(t_filtered,E_filtered,'color','black')
    title('Comparison between the compensated pulses and the fourier limit')
    text(500,0.75,sprintf('Curves have the\nsame integral'),'FontSize',5,'BackgroundColor',[.7 .9 .7])
    axis([-2000 2000 0 1])
    legend('fourier limit','smooth phase','super smooth','super smooth filtered')
    hold 'off'
    %---------------------------------subplot(2,2,2)-----------------------
    subplot(2,2,2)
    [ax,h1,h2] = plotyy([Ek_trebino(:,1) Ek(:,1)],[Ek_trebino(:,2) Ek(:,2)],[Ek_trebino(:,1) Ek(:,1)],[Ek_trebino(:,3) Ek(:,3)]);
    %legend([h1;h2],'intensity trebino','intensity modified','phase trebino', 'phase modified','Location','SouthOutside')
    xlabel('time in fs')
    title('Comparison between original pulse in time and the modified one')
    set(ax(1),'xlim',[-3000,4000])
    set(ax(2),'xlim',[-3000,4000])
    %---------------------------------subplot(2,2,3)-----------------------
    subplot(2,2,3)
    plot(cropped_original_omega.*1e-15,cropped_original_phase,'color','red')
    hold 'on'
    plot(cropped_original_omega.*1e-15,cropped_new_phase)
    plot(cropped_original_omega.*1e-15,(cropped_original_intensity .* 28)-15,'color','green')
    plot(omega,smoothed_GDD_integrated,'color','magenta')
    hold 'off'
    title('View of the whole phase curve and the smoothened area')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel('phase value')
    axis([1.80, 1.85, -15, 15])
    %---------------------------------subplot(2,2,4)-----------------------
    subplot(2,2,4)
    plot(omega,GDD_original_phase.*1e30, 'r')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    xlim([2*pi*c/1040e6, 2*pi*c/1018e6])
    hold 'on'
    plot(omega,GDD_new_phase.*1e30)
    plot(omega,smoothed.*1e30,'color','green')
    hold 'off'
    %-------------------------Figure(2)------------------------------------
    %----------------------------------------------------------------------
    figure(2)
    plot(omega,GDD_original_phase.*1e30, 'r')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    title('Comparison of all three GDD curves')
    axis([1.798 1.852 -1.5e6 1.5e6])
    hold 'on'
    plot(omega,GDD_new_phase.*1e30)
    plot(omega,smoothed.*1e30,'color','green')
    plot(omega,super_smoothed*1e30,'color','magenta','LineWidth',3)
    plot(omega,filtered_original*1e30,'color','black','LineWidth',3)
    hold 'off'
    legend('original','from smoothed phase','smoothed GDD from smoothed phase','super smooth')
    %-------------------------Figure(3)------------------------------------
    %----------------------------------------------------------------------
    figure(3)
    plot(cropped_original_omega.*1e-15,cropped_new_phase)
    hold on
    plot(omega,smoothed_GDD_integrated,'color',[0.5,0.5,0.5],'LineWidth',4)
    plot(omega,super_smoothed_integrated,'color','red')
    plot(original_omega*1e-15,original_phase,'color','cyan')
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_data_integrated_patched,'color','magenta')
    plot(cropped_original_omega.*1e-15,(cropped_original_intensity .* 22)-13,'color','green')
    hold off
    axis([1.80 1.85 -13 11])
    legend('smooth','smoothed GDD','super smooth','original','filtered GDD','spectrum')
    title('Direct comparison between the original phase curve and the integrated smoothed GDD (integration constants adapted)')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel('Phase')
    %-------------------------Figure(4)------------------------------------
    %----------------------------------------------------------------------
    figure(4)
    plot(omega,smoothed*1e30,'color','blue')
    axis([1.798 1.852 -1.2e6 1e6])
    title('Plain view of the super smoothed GDDs')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    hold 'on'
    plot(omega,super_smoothed*1e30,'color','green')
    hold 'off'
    legend('smooth','super smooth')
    %-------------------------Figure(5)------------------------------------
    %----------------------------------------------------------------------
    figure(5)
    plot(omega,super_smoothed*1e30,'color','green','LineWidth',3)
    hold on
    plot(omega,filtered_data*1e30,'r')
    plot(omega,filtered_original*1e30,'color','black')
    %plot(omega,temp*1e30,'color','magenta','LineWidth',3)
    hold off
    legend('super smoothed','super smoothed filtered','original filtered','zero phase filtering')
    title('Comparison of the best options for a simple GDD')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel(sprintf('GDD in fs^2'))
    xlim([1.798 1.85])
    %-------------------------Figure(6)------------------------------------
    %----------------------------------------------------------------------
    figure(6)
    plot(cropped_original_omega*1e-15,cropped_original_phase,'color','cyan')
    hold on
    plot(cropped_original_omega*1e-15,cropped_truncated_filtered_data_integrated_patched,'color','magenta')
    plot(cropped_original_omega*1e-15,cropped_diff_phase_filtered,'color','blue')
    plot(cropped_original_omega.*1e-15,(cropped_original_intensity .* 22)-13,'color','green')
    % Show the zero line
    plot(cropped_original_omega*1e-15,zeros(size(cropped_diff_phase_filtered)),'color','black')
    hold off
    axis([1.798 1.85 -13 11])
    legend('original','super smoothed filtered','difference','spectrum')
    title('Direct comparison between the original phase curve and the integrated smoothed GDD (integration constants adapted)')
    xlabel(sprintf('angular frequency in 10^{15} s^{-1}'))
    ylabel('Phase')
%%
    %saveas(gcf,sprintf('Daten/Pictures/sum_up/comparison_to_fourier_%s',endings{i}),'jpg')
end