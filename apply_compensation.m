clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
warning off
%format compact;
set(0,'DefaultFigureWindowStyle','docked')

folder = '../Daten/FROGS 10kHz_2mJ';
folderContent = dir(folder);
numberOfFiles = length(folderContent);
% Filtering of unwanted files necessary!
for i=1:numberOfFiles
    fprintf('File %d: %s\n',i,folderContent(i).name)
end

% 
% fprintf('---------RUNNING TEST SUITE---------\n')
% 
% fprintf('\nWe start with the following values:\n')
% fprintf('FOD: %2.2e\n',FOD*24)
% fprintf('TOD: %2.2e\n',TOD*6)
% fprintf('GDD: %2.2e\n',GDD*2)
% fprintf('GD: %2.2e\n',GD)
% fprintf('C: %2.2e\n',C)
% 
% options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
% min_func = @(x)compensation_minFuncForBruteForce_no_global2([FOD x(1) x(2) GD C],w0,w_Sk,I_Sk,p_Sk,Int_F);
% 
% [solution,val] = fminsearch(min_func,[TOD GDD],options);
% 
% fprintf('With the new optimization we find a value of %2.2f%% of the fourier peak\n',(1-val)*100)
% fprintf('The found values for TOD and GDD are:\n')
% fprintf('TOD: %2.2e\n',solution(1)*6)
% fprintf('GDD: %2.2e\n',solution(2)*2)
% 
% phase_new = polyval([FOD solution(1) solution(2) GD C],w_Sk-w0);
% N = 5*length(w_Sk);
% [w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
% [Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I_ext,l_ext);
% [Int_new,t_new,E_new] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
% factor = Int_F/Int_new;