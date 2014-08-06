%% ########################################################################
%
%  
%  Version: 3.0
%
%  INPUT: No input. At the moment the file path is hardcoded.
%
%  OUTPUT: Several plots showing the results of compression with varying 
%          degree.
%
%  HISTORY:
%      v1.0: First runable state (16.05.2014)
%      v2.0: Adjusted the code to give a comprehensive summary for the
%            mirror design. Removed most of the approximation schemes and
%            focus more on lower Taylor orders. (10.07.2014)
%      v2.1: Moved from a least square optimization method to a custom
%            method. This accommodates our needs better. (14.07.2014)
%      v2.2: Changed the method to only consider our spectral realm of
%            interest in order to avoid noise contribution. (15.07.2014)
%      v3.0: Concluded a series of test and code changes and ended up with
%            this running and stable version.
%
%% ########################################################################

%% INITILIZATION

    clc;    % Clear the command window.
    close all;  % Close all figures (except those of imtool.)
    clear;  % Erase all existing variables. Or clearvars if you want.
    format shortg;
    warning off
    %format compact;
    set(0,'DefaultFigureWindowStyle','docked')

    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
        cd(fileparts(which(mfilename)));
    end
    do_plot = false;
    max_order = 10;
    orders = 2:max_order;
%     orders = 4;
    
%%


%% Too make function definitions easier we define a bunch of global variables

    global I_Sk l_Sk p_Sk w_Sk
    global Int_F
    global p order figNum
    figNum = 1;

%%


%% This is the part where all filenames and directories are defined

    parent         = '../Daten/Chirped Mirrors/Frogs/25.3W_RTT=6.404us_Ip=12.8A/';
    filebase       = '001_AIR_FROG_25.3W_RTT=6.404us_Ip=12.8A.bin';  
    filename_Et    = sprintf('%s%s.Ek.dat',parent,filebase);
    filename_Speck = sprintf('%s%s.Speck.dat',parent,filebase);

%%


%% Open the file and get the original data
    
    % Time based field   
    Et   = dlmread(filename_Et);
    t_Et = Et(:,1);
    I_Et = Et(:,2);
    p_Et = Et(:,3);
    
    % Wavelength based field
    Sk   = dlmread(filename_Speck);
    l_Sk = Sk(:,1);
    I_Sk = Sk(:,2);
    p_Sk = Sk(:,3);

%%


%% Calculate the Fourier limit for our pulse

    Sk_cplx    = sqrt(I_Sk) .* exp(1i*0); % This equals constant phase -> Fourier limit
    [t_F,Ek_F] = Speck_Fourier(l_Sk*1e-9,Sk_cplx);
    t_F        = t_F * 1e15;
    Int_F      = abs(trapz(t_F,abs(Ek_F).^2));

%%


%% Physical units, transformations and physical calculations

    c    = 299792458;
    l0   = 1.030;
    w0   = 2*pi*0.299792458 / l0;
    w_Sk = 2*pi*c ./ (l_Sk*1e-9);
    w_Sk = w_Sk * 1e-15;
    
    % For differentiation it is good to have the x-values equally spaced
    w_Sk2 = linspace(w_Sk(1),w_Sk(end),length(w_Sk));
    I_Sk = interp1(w_Sk,I_Sk,w_Sk2);
    p_Sk = interp1(w_Sk,p_Sk,w_Sk2);
    w_Sk = w_Sk2; clear w_Sk2
    l_Sk = 2*pi*c ./ (w_Sk*1e15);
    l_Sk = l_Sk' * 1e9;
  
%%


%% Filter part

    % It turns out that it is beneficial to apply a low pass filter every
    % other differentiation in order to avoid the chaotic behaviour of
    % finite difference differentiation.
    sampling_freq   = 1/abs(w_Sk(2)-w_Sk(1));
    order           = 2;
    cut_off_freq    = 80;
    peak_to_peak_dB = 0.5;      % This is the Matlab doc recommend first guess
    stopband_atten  = 20;       % This is the Matlab doc recommend first guess
                                % The frequencies are always expressed normalized to the Nyquist 
                                % frequency, i.e. half the sampling rate
    normed_cutoff   = cut_off_freq / (0.5*sampling_freq);
                                % The following function call creates the filter fractions given in the
                                % form [nominators,denominators]
    [B,A]           = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
    %-------------------------filter the data------------------------------
                                % filtfilt has the big advantage over filter that it operates at a zero phase
                                % shift, thus allowing to obtain a direct correspondance to the filtered data
    filtered_p_Sk   = filtfilt(B,A,p_Sk);
    
    % Use this to find our spectral area of interest
    M        = max(I_Sk);
    perc     = 0.05;
    fit_w_Sk = w_Sk(I_Sk > perc*M);
    % Now find the fringe values
    lower  = min(fit_w_Sk);
    higher = max(fit_w_Sk);
    % Now we use these values to include everything in between
    fit_w_Sk = w_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_p_Sk = filtered_p_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_I_Sk = I_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_l_Sk = l_Sk(w_Sk >= lower & w_Sk <= higher);
    
    % Make a second Fourier limit to proof that our chosen spectral realm
    % of interest. The integrals should be the same.

    Sk_cplx   = sqrt(fit_I_Sk) .* exp(1i*0); % This equals constant phase -> Fourier limit
    lambda    = 2*pi*c ./ (fit_w_Sk*1e15);
    [fit_t_F,fit_Ek_F] = Speck_Fourier(lambda,Sk_cplx);
    fit_t_F   = fit_t_F * 1e15;
    fit_Int_F = abs(trapz(t_F,abs(Ek_F).^2));
    
%%


%% Taylor part

    D1 = diff(filtered_p_Sk) ./ (w_Sk(1)-w_Sk(2));
    % Immediately apply our lowpass filter
    %D1 = filtfilt(B,A,D1);
    D2 = diff(D1) ./ (w_Sk(1)-w_Sk(2));
    D2 = filtfilt(B,A,D2);
    D3 = diff(D2) ./ (w_Sk(1)-w_Sk(2));
    %D3 = filtfilt(B,A,D3);
    
    % Get our Taylor coefficients T1,T2,...
    T0 = interp1(w_Sk,filtered_p_Sk,w0);
    T1 = interp1(w_Sk(1:end-1),D1,w0);
    T2 = interp1(w_Sk(1:end-2),D2,w0);
    T3 = interp1(w_Sk(1:end-3),D3,w0);
    
    % Get our Taylor approximation of order O(0),O(1),...
    % This proves to bring no plausible results since we merely
    % differentiate the numerical data. In this case the concept of the
    % Taylor series is broken! We don't have an analytical function...
    O0 = T0*ones(size(w_Sk));
    O1 = O0 + T1*(w_Sk-w0);
    O2 = O1 + T2/2*(w_Sk-w0).^2;
    O3 = O2 + T3/6*(w_Sk-w0).^3;
%%


%% Fit Polynomials to the phase and plot the results depending on degree

    % These vectors will hold our results while we step through the orders
    Pv = zeros(max_order,length(w_Sk));
    P_optv = zeros(max_order,length(w_Sk));
    D2_optv = zeros(max_order,length(w_Sk)-2);
    D3_optv = zeros(max_order,length(w_Sk)-3);
    
    Ev = zeros(max_order,length(w_Sk));
    E_optv = zeros(max_order,length(w_Sk));
    opt_peaks = zeros(1,max_order);
    opt_best_order = 0;
    opt_best_value = 0;
    peaks = zeros(1,max_order);
    best_order = 0;
    best_value = 0;
    for order=orders
        p = polyfit(fit_w_Sk,fit_p_Sk,order);
        P = polyval(p,w_Sk);

        % Matlab's builtin `polyfit` uses least square optimization. This might
        % be an inadequate measure for our Fourier transform needs.
        % Thus we have to write our own minimization function and use Matlab's
        % built in `fminsearch`.
        [solution,val] = make_fourier_fit();
        P_opt = polyval(solution,w_Sk);
        D1_opt = diff(P_opt) ./ (w_Sk(1)-w_Sk(2));
        D2_opt = diff(D1_opt) ./ (w_Sk(1)-w_Sk(2));
        D2_opt = filtfilt(B,A,D2_opt);
        D3_opt = diff(D2_opt) ./ (w_Sk(1)-w_Sk(2));
        Pv(order,:) = P;
        P_optv(order,:) = P_opt;
        D2_optv(order,:) = D2_opt;
        D3_optv(order,:) = D3_opt;

%%


%% Integrate back into the time domain

        d     = p_Sk-Pv(order,:);
        d_opt = p_Sk-P_optv(order,:);

        Sk_cplx     = sqrt(I_Sk) .* exp(1i*d);
        Sk_cplx_opt = sqrt(I_Sk) .* exp(1i*d_opt);

        [t,E]         = Speck_Fourier(l_Sk'*1e-9,Sk_cplx);
        [t_opt,E_opt] = Speck_Fourier(l_Sk'*1e-9,Sk_cplx_opt);

        t     = t * 1e15;
        t_opt = t_opt * 1e15;

        Int     = abs(trapz(t,abs(E).^2));
        Int_opt = abs(trapz(t_opt,abs(E_opt).^2));

        % Scale the integrals to allow for comparison with Fourier limit
        E     = E .* sqrt(Int_F ./ Int);
        E_opt = E_opt .* sqrt(Int_F ./ Int_opt);
        fprintf('Order %d: The optimization for the whole range leads to %2.2f percent of the Fourier limit.\n',order,max(abs(E_opt).^2)*100)
        Ev(order,:) = E;
        E_optv(order,:) = E_opt;
        opt_peaks(1,order) = max(abs(E_opt).^2)*100;
        if max(abs(E_opt).^2)*100 > opt_best_value
            opt_best_order = order;
            opt_best_value = max(abs(E_opt).^2)*100;
        end
        peaks(1,order) = max(abs(E).^2)*100;
        if max(abs(E).^2)*100 > best_value
            best_order = order;
            best_value = max(abs(E).^2)*100;
        end
        
    figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(t_F,abs(Ek_F).^2,'b')
        plot(t,abs(E).^2,'g')
        % Shift it correctly
        t0 = find_closest_idx(t_opt,0);
        E0 = find_closest_idx(abs(E_opt).^2,max(abs(E_opt).^2));
        plot(t_opt,circshift(abs(E_opt).^2,t0-E0),'c')
        hold off
        title('Achieved compression compared to Fourier limit')
        xlabel('Time[fs]')
        xlim([-2000 2000])
        legend('Fourier limit','Least square','Optimization')
    figure(figNum)
            figNum = figNum + 1;
            hold on
            plot(w_Sk,I_Sk./max(I_Sk).*max(filtered_p_Sk),'k')
            plot(w_Sk,p_Sk,'b')
            plot(w_Sk,P,'r')
            plot(w_Sk,P_opt,'g')
            plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
            hold off
            xlim([lower-0.01 higher+0.01])
            ylim([min(p_Sk)-2 max(p_Sk)+2])
            title('Taylor approximations for our phase curve')
            xlabel('Omega[1/fs]')
            legend('Intensity','Original Phase','least square','custom')
    end
%%


%% Figure part
    if do_plot
        %-------------------------Figure(1)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        plotyy(t_Et,I_Et,t_Et,p_Et)
        title('Intensity and phase of our pulse in the time domain')
        xlabel('Time[fs]')
        legend('Intensity','Phase')
        %-------------------------Figure(2)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        plotyy(w_Sk,I_Sk,w_Sk,p_Sk)
        title('Intensity and phase of our pulse in the frequency domain')
        xlabel('Omega[1/fs]')
        legend('Intensity','Phase')
        %-------------------------Figure(3)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk,p_Sk)
        plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
        hold off
        title('Phase in spectral domain')
        xlabel('Omega[1/fs]')
        %-------------------------Figure(4)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk(1:end-1),D1)
        plot([w0 w0],[min(D1) max(D1)],'r')
        hold off
        title('GD of our pulse (filter applied)')
        xlabel('Omega[1/fs]')
        %-------------------------Figure(5)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk(1:end-2),D2)
        plot(w_Sk(1:end-2),D2_optv(opt_best_order,:),'g')
        plot([w0 w0],[min(D2) max(D2)],'r')
        hold off
        title('GDD of our pulse (filter applied)')
        xlabel('Omega[1/fs]')
        legend('filtered phase','custom')
        %-------------------------Figure(6)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk(1:end-3),D3)
        plot(w_Sk(1:end-3),D3_optv(opt_best_order,:),'g')
        plot([w0 w0],[min(D3) max(D3)],'r')
        hold off
        title('TOD of our pulse (filter applied)')
        xlabel('Omega[1/fs]')
        legend('filtered phase','custom')
        %-------------------------Figure(7)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk,I_Sk./max(I_Sk).*max(filtered_p_Sk),'k')
        plot(w_Sk,p_Sk,'b')
        plot(w_Sk,Pv(opt_best_order,:),'r')
        plot(w_Sk,P_optv(opt_best_order,:),'g')
        plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
        hold off
        xlim([lower-0.01 higher+0.01])
        ylim([min(p_Sk)-2 max(p_Sk)+2])
        title('Taylor approximations for our phase curve')
        xlabel('Omega[1/fs]')
        legend('Intensity','Original Phase','least square','custom')
        %-------------------------Figure(8)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(w_Sk,I_Sk./max(I_Sk).*50,'k')
        plot(w_Sk,p_Sk-Pv(opt_best_order,:),'g')
        plot(w_Sk,p_Sk-P_optv(opt_best_order,:),'c')
        plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
        hold off
        xlim([lower-0.01 higher+0.01])
        ylim([-20 50])
        title('Difference between the Taylor approximation and real phase')
        xlabel('Omega[1/fs]')
        legend('Intensity','Original Phase','least square','custom')
        %-------------------------Figure(9)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(t_F,abs(Ek_F).^2,'b')
        plot(t,abs(E).^2,'g')
        % Shift it correctly
        t0 = find_closest_idx(t_opt,0);
        E0 = find_closest_idx(abs(E_optv(opt_best_order,:)).^2,max(abs(E_optv(opt_best_order,:)).^2));
        plot(t_opt,circshift(abs(E_optv(opt_best_order,:)).^2,t0-E0),'c')
        hold off
        title('Achieved compression compared to Fourier limit')
        xlabel('Time[fs]')
        xlim([-2000 2000])
        legend('Fourier limit','Least square','Optimization')
        %-------------------------Figure(10)------------------------------------
        %----------------------------------------------------------------------
        figure(figNum)
        figNum = figNum + 1;
        hold on
        plot(1:max_order,opt_peaks)
        plot(1:max_order,peaks,'r')
        hold off
        title('Compression quality with increasing order')
        xlabel('Order')
        ylabel('Peak percentage of Fourier limit')
        legend('custom','least squares')
    end
%%


%% Clear all unwanted variables

    clear Sk
    clear vq
    %clear increment start ende
%%
    %saveas(gcf,sprintf('Daten/Pictures/sum_up/comparison_to_fourier_%s',endings{i}),'jpg')