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
%      v1.0: First runable state (12.08.2014)
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

global I_Sk l_Sk p_Sk w_Sk
global Int_F
global figNum
figNum = 1;

%% Start Algorithm

[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData();

[Int_F,t_F,Ek_F]                      = compensation_calcFourierlimit(I_Sk,l_Sk);

[l0,w0,w_Sk,I_Sk,p_Sk,l_Sk]           = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);

[B,A]                                 = compensation_makeFilter(w_Sk);

filtered_p_Sk                         = filtfilt(B,A,p_Sk);

[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,filtered_p_Sk);

[D1,D2,D3]                            = compensation_calcDevs(filtered_p_Sk,w_Sk,B,A);

% Brute force all the polynomials of the form ax^3+bx^2+cx+d, where c and d
% remain constant

% Find a start hint
p0 = polyfit(fit_w_Sk,fit_p_Sk,3);

% Prepare the range
orda = floor(log10(p0(1)));
ordb = floor(log10(p0(2)));

% Generate the coefficients
%p = ...

% Calculate the value
val = compensation_minFuncForBruteForce(p);

% Store values for later plotting
% a(*) = a;
% b(*) = b;
% v(*) = v;

% Make a 2D plot to visualize the result























