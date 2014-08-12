%% ########################################################################
%
%  
%  Version: 1.0
%
%  INPUT: No input. At the moment the file path is hardcoded.
%
%  OUTPUT: Several plots showing the results of compression with varying 
%          degree.
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
    warning off
    %format compact;
    set(0,'DefaultFigureWindowStyle','docked')

    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
        cd(fileparts(which(mfilename)));
    end
    do_plot = false;
    
%%


%% Too make function definitions easier we define a bunch of global variables

    global I_Sk l_Sk p_Sk w_Sk
    global Int_F
    global p p_bkp order figNum
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
    
    % Use this to find our spectral area of interest
    M        = max(I_Sk);
    perc     = 0.05;
    fit_w_Sk = w_Sk(I_Sk > perc*M);
    % Now find the fringe values
    lower  = min(fit_w_Sk);
    higher = max(fit_w_Sk);
    % Now we use these values to include everything in between
    fit_w_Sk = w_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_p_Sk = p_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_I_Sk = I_Sk(w_Sk >= lower & w_Sk <= higher);
    fit_l_Sk = l_Sk(w_Sk >= lower & w_Sk <= higher);
       
%%


%% Taylor part

%     D1 = diff(filtered_p_Sk) ./ (w_Sk(1)-w_Sk(2));
%     % Immediately apply our lowpass filter
%     %D1 = filtfilt(B,A,D1);
%     D2 = diff(D1) ./ (w_Sk(1)-w_Sk(2));
%     D2 = filtfilt(B,A,D2);
%     D3 = diff(D2) ./ (w_Sk(1)-w_Sk(2));
%     %D3 = filtfilt(B,A,D3);

%%


%% Fit Polynomials to the phase and plot the results depending on degree

    order = 4;
    p = polyfit(fit_w_Sk,fit_p_Sk,order);
    P = polyval(p,w_Sk);
    
    % Make our initial calculations as usual
    d       = p_Sk-P;
    Sk_cplx = sqrt(I_Sk) .* exp(1i*d);
    [t,E]   = Speck_Fourier(l_Sk'*1e-9,Sk_cplx);
    t       = t * 1e15;
    Int     = abs(trapz(t,abs(E).^2));
    % Scale the integrals to allow for comparison with Fourier limit
    E       = E .* sqrt(Int_F ./ Int);

%%


%% Here we change the starting polynomial

    fraction = 0.00001; % This is the limit. If we increase fraction further our algorithm converges.
    p_bkp = p;
    
figure(1)
    hold on
    plot(t_F,abs(Ek_F).^2,'b')
    plot(t,abs(E).^2,'g')
    hold off
    title('Achieved compression compared to Fourier limit')
    xlabel('Time[fs]')
    ylabel('Arbitrary units')
    xlim([-500 500])
    ylim([0 1])
    legend('Fourier limit','Least square','Optimization')
    
figure(2)
    hold on
    plot(w_Sk,I_Sk./max(I_Sk).*max(p_Sk),'k')
    plot(w_Sk,p_Sk,'b')
    plot(w_Sk,P,'r')
    hold off
    xlim([lower-0.01 higher+0.01])
    ylim([-30 50])
    title('Taylor approximations for our phase curve')
    xlabel('Omega[1/fs]')
    ylabel('Phase[rad]')
    legend('Intensity','Original Phase','least square','custom')
        
figure(3)
    plot(w_Sk, d)
    ylim([-50 50])
    title('Difference between the original phase curve and the found optima')
    xlabel('Omega[1/fs]')
    ylabel('Difference[rad]')

    
for i=1:10
    rands = (rand([1,order+1])-0.5)*10;
    rands = rands.*fraction;
    p = p_bkp + rands.*p_bkp;
    disp(p)
    
    % Matlab's builtin `polyfit` uses least square optimization. This might
    % be an inadequate measure for our Fourier transform needs.
    % Thus we have to write our own minimization function and use Matlab's
    % built in `fminsearch`.
    [solution,val] = compensation_improveFit();
    P_opt = polyval(solution,w_Sk);

    d_opt = p_Sk-P_opt;
    Sk_cplx_opt = sqrt(I_Sk) .* exp(1i*d_opt);
    [t_opt,E_opt] = Speck_Fourier(l_Sk'*1e-9,Sk_cplx_opt);
    t_opt = t_opt * 1e15;
    Int_opt = abs(trapz(t_opt,abs(E_opt).^2));
    % Scale the integrals to allow for comparison with Fourier limit
    E_opt = E_opt .* sqrt(Int_F ./ Int_opt);

figure(1)
    hold on
    % Shift it correctly
    t0 = find_closest_idx(t_opt,0);
    E0 = find_closest_idx(abs(E_opt).^2,max(abs(E_opt).^2));
    plot(t_opt,circshift(abs(E_opt).^2,[1,t0-E0]),'c')
    hold off
figure(2)
    hold on
    plot(w_Sk,P_opt,'g')
    hold off
figure(3)
    hold on
    plot(w_Sk,d_opt)
    hold off
end