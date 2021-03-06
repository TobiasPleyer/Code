clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
format longg;
warning off
%format compact;
set(0,'DefaultFigureWindowStyle','docked')

pattern = '^\d\d\d_';

%% Simple pulse
folder = '../Daten/FROGS 10kHz_2mJ/New/';
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
    error('No files found that match the patter!')
end

factors = zeros(size(usefulFiles));
figure()
hold on
for i=1:length(usefulFiles)
    fileBase = usefulFiles{i};
    [t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData(folder,fileBase);
    [l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk,1030);
    [fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);
    temp_polynomial                       = polyfit(fit_w_Sk,fit_p_Sk,4);
    GDD = polyval(polyder(polyder(temp_polynomial)),w0) / 2;
    if GDD > 0
        factors(i) = 1;
    else
        factors(i) = -1;
    end
    binary = dec2bin(i);
    switch length(binary)
        case 1
            binary = ['00' binary];
        case 2
            binary = ['0' binary];
    end
    plot(fit_w_Sk,factors(i)*fit_p_Sk,'Color',[str2num(binary(1)) str2num(binary(2)) str2num(binary(3))])
end
xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
ylabel('[rad]','fontweight','bold','fontsize',16);
title('Direct comparison between the phase curves for our 6 retrievals','fontweight','bold','fontsize',16);
hold off

fprintf('\n--------------RUNNING TEST SUITE--------------\n')

% We take the first file in our list, find its values for GD, GDD, TOD,
% etc.. and find an optimal solution. When this is done we load the other
% files and see how well we can compress them with that solution.
fprintf('Opening first file...\n')
fileBase = usefulFiles{6};
[w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,factors(6),1030,'');
    
fprintf('\nApplying the found optimum to the remaining pulses\n')

for i=1:length(usefulFiles)
    fileBase = usefulFiles{i};
    [w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,factors(i),1030,sprintf('for the file:\n%s',fileBase),phase_new);
end
c    = 299.792458;
l_Sk = 2*pi*c./w_Sk;
I_Sk = I_Sk./l_Sk.^2;
I_Sk = I_Sk./max(I_Sk);

l_export = linspace(l_Sk(1),l_Sk(end),1000);
p_export = interp1(l_Sk,phase_new,l_export,'spline');
I_export = interp1(l_Sk,I_Sk,l_export);

% dlmwrite('../Daten/FROGS 10kHz_2mJ/Trubetskov_simple_Pulse_Approximation_Phase.dat',[l_export',p_export'],'precision',8)
% dlmwrite('../Daten/FROGS 10kHz_2mJ/Trubetskov_simple_Pulse_Approximation_Intensity.dat',[l_export',I_export'],'precision',8)


%% Complex pulse
folder = '../Daten/FROGS 10kHz_2mJ/4mJ/';
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
    error('No files found that match the patter!')
end
for i=1:length(usefulFiles)
    fileBase = usefulFiles{i};
    [w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,1,1030,sprintf('for the file:\n%s',fileBase));
end

% dlmwrite('../Daten/FROGS 10kHz_2mJ/Trubetskov_complex_Pulse_original_Phase.dat',[w_Sk',p_Sk',I_Sk'])