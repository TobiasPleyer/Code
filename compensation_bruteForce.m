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

[l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);

[B,A]                                 = compensation_makeFilter(w_Sk);

filtered_p_Sk                         = filtfilt(B,A,p_Sk);

[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);

[D1,D2,D3]                            = compensation_calcDevs(filtered_p_Sk,w_Sk,B,A);

% Brute force all the polynomials of the form ax^3+bx^2+cx+d, where c and d
% remain constant

% Find a start hint
p0 = polyfit(fit_w_Sk,fit_p_Sk,3);

% Prepare the range
ord_a = floor(log10(abs(p0(1))));
ord_b = floor(log10(abs(p0(2))));

% Generate the coefficients
n_a    = 1*10^(ord_b+3);
n_a    = 10000;
min_a  = 1*10^(ord_a-3);
max_a  = 9*10^(ord_a+1);
step_a = (max_a - min_a)/n_a;

n_b    = 1*10^(ord_b+3);
n_b    = 10000;
min_b  = 1*10^(ord_b-3);
max_b  = 9*10^(ord_b+1);
step_b = (max_b - max_b)/n_b;

% Reserve arrays for storage
A = linspace(min_a,max_a,n_a);
B = linspace(min_b,max_b,n_b);
V = zeros(n_b,n_a);

% Enter a loop to brute force the solutions
N = 1;
for a=A
    for b=B
        p = [a b 0 0];
        % Calculate the value
        val = compensation_minFuncForBruteForce(p);
        V(N) = val;
        if mod(N,1000)==0
            fprintf('N: %d\n',N)
        end
        N = N + 1;
    end
end

% For longer runs to save the data for later revision
save('Variables/run03')
% figure(figNum)
% figNum = figNum + 1;
% % Make a 3D surf plot to visualize the result
% h = surf(A,B,V);
% % Necessary, or the lines of the grid will overlay the colors
% set(h,'LineStyle','none')






















