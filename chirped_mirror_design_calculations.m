%% Clear and prepare the workspace
    clear
    format shortg
    set(0,'DefaultFigureWindowStyle','docked')

%% Define common variables
    mirror    = xlsread('T:/Tobias/Daten/Chirped Mirrors/First design/PC526_Result_4b_Spect_Data_GDD_Est.xlsx','Sheet2');
    original  = dlmread('T:/Tobias/Daten/Chirped Mirrors/First design/first_design_equidistant_data_of_pulse.csv');
    new_pulse = dlmread('T:/Tobias/Daten/Chirped Mirrors/Frogs/AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A/001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A.bin.Speck.dat','\t');
    c         = 299792458;

    %----Michael Trubetskov's curves-------------------------------------------
    m_wavelength = mirror(:,1);
    m_phase      = mirror(:,3);
    m_GD         = mirror(:,4);
    m_GDD        = mirror(:,5);
    m_omega      = mirror(:,6);
    m_GDD_num    = mirror(:,7);

    %----Original data curves--------------------------------------------------
    o_wavelength = original(:,1);
    o_omega      = 2*pi*c ./ (o_wavelength*1e-9);
    o_phase      = original(:,3);
    o_phase = o_phase((o_omega > 1.8e15) & (o_omega < 1.85e15));
    o_omega = o_omega((o_omega > 1.8e15) & (o_omega < 1.85e15));

    %----Newest pulse data-----------------------------------------------------
    p_wavelength = new_pulse(:,1);
    p_omega      = 2*pi*c ./ (p_wavelength*1e-9);
    p_omega      = linspace(p_omega(1),p_omega(end),length(p_omega));
    p_phase      = new_pulse(:,3);
    p_phase = p_phase((p_omega > 1.8e15) & (p_omega < 1.85e15));
    p_omega = p_omega((p_omega > 1.8e15) & (p_omega < 1.85e15));

    %----Final touches---------------------------------------------------------
    s = smooth(o_phase);
    o_omega = o_omega * 1e-15;
    p_omega = p_omega * 1e-15;
    %o_phase = o_phase + 720*o_omega - 1300;
    clear num new_pulse mirror p_omega2 original o_omega2;

%% Low pass filtering
    %--------Filter the phase data-----------------------------------------
    %-------------------------Make preparations for a lowpass filter-------
    sampling_freq = (-1)/(o_omega(2)-o_omega(1)); % 0.2THz here
    order         = 2;
    cut_off_freq  = 80;
    peak_to_peak_dB = 0.5;      % This is the Matlab doc recommend first guess
    stopband_atten  = 20;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A] = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
    %-------------------------filter the data------------------------------
                                % filtfilt has the big advantage over filter that it operates at a zero phase
                                % shift, thus allowing to obtain a direct correspondance to the filtered data
    o_filtered_phase  = filtfilt(B,A,o_phase);
    p_filtered_phase  = filtfilt(B,A,p_phase);
    
%% Calculation part

    %----Derive from unfiltered phase--------------------------------------
    o_GD = diff(o_phase)./(o_omega(1)-o_omega(2));
    o_GDD = diff(smooth(o_GD))./(o_omega(2)-o_omega(1));
    %----Derive from original filtered phase-------------------------------
    o_filtered_GD  = diff(o_filtered_phase)./(o_omega(1)-o_omega(2));
    o_filtered_GDD = diff(o_filtered_GD)./(o_omega(2)-o_omega(1));
    %--------A final filter is advisable-----------------------------------
    cut_off_freq = 250; normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
    [B,A] = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
    o_filtered_GDD2  = filtfilt(B,A,o_filtered_GDD);
    %----Derive from unfiltered phase--------------------------------------
    %----------------------------------------------------------------------
    p_GD = diff(p_phase)./(p_omega(1)-p_omega(2));
    p_GDD = diff(smooth(p_GD))./(p_omega(2)-p_omega(1));
    %----Derive from new filtered phase-------------------------------
    p_filtered_GD  = diff(p_filtered_phase)./(p_omega(1)-p_omega(2));
    p_filtered_GDD = diff(p_filtered_GD)./(p_omega(2)-p_omega(1));
    %--------A final filter is advisable-----------------------------------
    p_filtered_GDD2  = filtfilt(B,A,p_filtered_GDD);
    
    %----Find the difference between the reference and the new pulse-------
    vq     = interp1(o_omega,o_phase,p_omega);
    D1     = p_phase - vq';
    vq     = interp1(o_omega,o_filtered_phase,p_omega);
    D2     = p_filtered_phase - vq';
    GD_D2  = diff(D2) ./ (p_omega(1)-p_omega(2));
    GDD_D2 = diff(GD_D2) ./ (p_omega(1)-p_omega(2));
%% Plotting part

    %-------------------------Figure(1)------------------------------------
    %----------------------------------------------------------------------
    figure(1)
    plot(m_omega,m_phase)
    hold on
    plot(o_omega,o_phase,'r')
    plot(p_omega,-p_phase,'m')
    plot(o_omega,o_filtered_phase,'g')
    plot(p_omega,-p_filtered_phase,'k')
    hold off
    title('Phase curve from Michael Trubetskov')
    legend('Trubetskov','Reference Pulse','New Pulse')
    xlabel('omega [fs^-1]')
    ylabel('phase [rad]')
    %-------------------------Figure(2)------------------------------------
    %----------------------------------------------------------------------
    figure(2)
    plot(m_omega,m_GD)
    hold on
    plot(o_omega(1:end-1),-o_filtered_GD,'r')
    plot(p_omega(1:end-1),-p_filtered_GD,'g')
    hold off
    title('GD curves')
    legend('Trubetskov','reference','new')
    xlabel('omega [fs^-1]')
    ylabel('GD [fs]')
    %-------------------------Figure(3)------------------------------------
    %----------------------------------------------------------------------
    figure(3)
    plot(m_omega,m_GDD)
    hold on
    plot(o_omega(1:end-2),-o_filtered_GDD2,'r')
    plot(p_omega(1:end-2),-p_filtered_GDD2,'g')
    hold off
    title('GDD curves')
    legend('Trubetskov','old','new')
    xlabel('omega [1/fs]')
    ylabel('GDD [fs^2]')
    %-------------------------Figure(4)------------------------------------
    %----------------------------------------------------------------------
    figure(4)
    plot(p_omega,D1)
    hold on
    plot(p_omega,D2,'r')
    hold off
    title('Differnce between the reference pulse and the new pulse')
    legend('raw','filtered')
    xlabel('omega [1/fs]')
    ylabel('phase [rad]')
    %-------------------------Figure(5)------------------------------------
    %----------------------------------------------------------------------
    figure(5)
    plot(p_omega(1:end-1),GD_D2)
    title('GD of the difference old-new')
    xlabel('omega [1/fs]')
    %-------------------------Figure(6)------------------------------------
    %----------------------------------------------------------------------
    figure(6)
    plot(p_omega(1:end-2),filtfilt(B,A,GDD_D2))
    title('GDD of the difference old-new')
    xlabel('omega [1/fs]')