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
%      v4.0: Modularized the main code to make heavy use of small
%            specialized function calls. This ensures readability,
%            flexibility and maintainability. (12.08.2014)
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

global I_Sk l_Sk p_Sk w_Sk
global Int_F
global p order figNum
figNum = 1;

%% Start Algorithm

[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData();

[Int_F,t_F,Ek_F]                      = compensation_calcFourierlimit(I_Sk,l_Sk);

[l0,w0,w_Sk,I_Sk,p_Sk,l_Sk]           = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);

[B,A]                                 = compensation_makeFilter(w_Sk);

filtered_p_Sk                         = filtfilt(B,A,p_Sk);

[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,filtered_p_Sk);

[D1,D2,D3]                            = compensation_calcDevs(filtered_p_Sk,w_Sk,B,A);

%% Fit Polynomials to the phase and plot the results depending on degree

% These vectors will hold our results while we step through the orders
Pv             = zeros(max_order,length(w_Sk));
P_optv         = zeros(max_order,length(w_Sk));
D2_optv        = zeros(max_order,length(w_Sk)-2);
D3_optv        = zeros(max_order,length(w_Sk)-3);
Ev             = zeros(max_order,length(w_Sk));
E_optv         = zeros(max_order,length(w_Sk));
opt_peaks      = zeros(1,max_order);
opt_best_order = 0;
opt_best_value = 0;
peaks          = zeros(1,max_order);
best_order     = 0;
best_value     = 0;

for order=orders
    p = polyfit(fit_w_Sk,fit_p_Sk,order);
    P = polyval(p,w_Sk);

    % Matlab's builtin `polyfit` uses least square optimization. This might
    % be an inadequate measure for our Fourier transform needs.
    % Thus we have to write our own minimization function and use Matlab's
    % built in `fminsearch`.
    [solution,val]         = compensation_improveFit();
    P_opt                  = polyval(solution,w_Sk);
    [D1_opt,D2_opt,D3_opt] = compensation_calcDevs(P_opt,w_Sk,B,A);
    Pv(order,:)            = P;
    P_optv(order,:)        = P_opt;
    D2_optv(order,:)       = D2_opt;
    D3_optv(order,:)       = D3_opt;

    [t,E]                  = compensation_calcCompression(p_Sk,Pv(order,:),I_Sk,l_Sk,Int_F);
    [t_opt,E_opt]          = compensation_calcCompression(p_Sk,P_optv(order,:),I_Sk,l_Sk,Int_F);
    
    fprintf('Order %d: The optimization for the whole range leads to %2.2f percent of the Fourier limit.\n',order,max(abs(E_opt).^2)*100)
    Ev(order,:)            = E;
    E_optv(order,:)        = E_opt;
    opt_peaks(1,order)     = max(abs(E_opt).^2)*100;
    if max(abs(E_opt).^2)*100 > opt_best_value
        opt_best_order = order;
        opt_best_value = max(abs(E_opt).^2)*100;
    end
    peaks(1,order)         = max(abs(E).^2)*100;
    if max(abs(E).^2)*100 > best_value
        best_order = order;
        best_value = max(abs(E).^2)*100;
    end

    compensation_plot2(w0    ...
                      ,t_F   ...
                      ,Ek_F  ...
                      ,t_opt ...
                      ,E_opt ...
                      ,t     ...
                      ,E     ...
                      ,w_Sk  ...
                      ,I_Sk  ...
                      ,p_Sk  ...
                      ,filtered_p_Sk ...
                      ,P     ...
                      ,P_opt)
end

if do_plot
    compensation_plot(w0 ...
                     ,t_Et ...
                     ,I_Et ...
                     ,p_Et ...
                     ,w_Sk ...
                     ,I_Sk ...
                     ,p_Sk ...
                     ,filtered_p_Sk ...
                     ,D1   ...
                     ,D2   ...
                     ,D3   ...
                     ,D2_optv(opt_best_order,:) ...
                     ,D3_optv(opt_best_order,:) ...
                     ,Pv(opt_best_order,:)      ...
                     ,P_optv(opt_best_order,:)  ...
                     ,t_F   ...
                     ,Ek_F  ...
                     ,t     ...
                     ,E     ...
                     ,t_opt ...
                     ,E_optv(opt_best_order,:) ...
                     ,max_order ...
                     ,opt_peaks ...
                     ,peaks)
end

clear Sk
clear vq