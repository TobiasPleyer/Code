clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
warning off
%format compact;
set(0,'DefaultFigureWindowStyle','docked')

pattern = '^\d\d_';
folder = '../Daten/FROGS 10kHz_2mJ/Old/';
folderContent = dir(folder);
numberOfFiles = length(folderContent);
usefulFiles = {};
% Filtering of unwanted files necessary!
fprintf('Searching for files with pattern: %s\n\n',pattern)
for i=1:numberOfFiles
    fileName = folderContent(i).name;
    if ~isempty(regexp(fileName,pattern,'match'))
        fileName = strrep(fileName,'.Speck.dat','');
        fileName = strrep(fileName,'.Ek.dat','');
        L = length(usefulFiles);
        flag = true;
        for j=1:L
            if fileName==usefulFiles{j}
                flag = false;
            end
        end
        if flag
            usefulFiles{L+1} = fileName;
            fprintf('File %d: %s\n',L+1,fileName)
        end
    end
end
if isempty(usefulFiles)
    fprintf('No files found that match the patter!\n')
end

figure(1)
hold on
for i=1:length(usefulFiles)
    fileBase = usefulFiles{i};
    [t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData(folder,fileBase);
    [l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);
    [fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);
    plot(fit_w_Sk,fit_p_Sk)
end
hold off

fprintf('\n--------------RUNNING TEST SUITE-------------\n')

% We take the first file in our list, find its values for GD, GDD, TOD,
% etc.. and find an optimal solution. When this is done we load the other
% files and see how well we can compress them with that solution.
fprintf('Opening first file...\n')
fileBase = usefulFiles{5};
[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData(folder,fileBase);
[Int_F,t_F,Ek_F]                      = compensation_calcFourierlimit(I_Sk,l_Sk);
[l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);
[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);
polynomial = polyfit(fit_w_Sk,fit_p_Sk,4);
% Find the correct coefficients
C = polyval(polynomial,w0);
GD = polyval(polyder(polynomial),w0);
GDD = polyval(polyder(polyder(polynomial)),w0) / 2;
TOD = polyval(polyder(polyder(polyder(polynomial))),w0) / 6;
FOD = polyval(polyder(polyder(polyder(polyder(polynomial)))),w0) / 24;

fprintf('\nWe start with the following values found from our least square approximation:\n')
fprintf('TOD: %2.2e\n',TOD*6)
fprintf('GDD: %2.2e\n',GDD*2)
fprintf('GD: %2.2e\n',GD)
fprintf('C: %2.2e\n',C)

options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
min_func = @(x)compensation_minFuncForBruteForce_no_global2([x(1) x(2) GD C],w0,w_Sk,I_Sk,p_Sk,Int_F);

[solution,val] = fminsearch(min_func,[TOD GDD],options);

fprintf('\nWith our optimization routine we find a value of %2.2f%% of the fourier peak\n',(1-val)*100)
fprintf('The found values for TOD and GDD are:\n')
fprintf('TOD: %2.2e\n',solution(1)*6)
fprintf('GDD: %2.2e\n',solution(2)*2)

phase_new = polyval([solution(1) solution(2) GD C],w_Sk-w0);
N = 5*length(w_Sk);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I_ext,l_ext);
[Int_new,t_new,E_new] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor = Int_F/Int_new;

figure(2)
    hold on
    plot(w_Sk,p_Sk,'r');
    plot(w_Sk,phase_new,'b');
    legend('Original','New found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 7: Comparison of the new found optimum to the original phase','fontweight','bold','fontsize',16);
    hold off
    
figure(3)
    plot(w_Sk,p_Sk-phase_new);
    legend('Original phase minus new found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 7: Difference between original phase and new optimum','fontweight','bold','fontsize',16);

figure(4)
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t_new,abs(E_new).^2.*factor,'r')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit','Compressed from new optimum')
    title('Test 7: Observation of the compression of our new found optimum','fontweight','bold','fontsize',16);