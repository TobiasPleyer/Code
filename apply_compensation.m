clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
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
    [l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);
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
[w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,factors(6));
    
fprintf('\nApplying the found optimum to the remaining pulses\n')

for i=1:length(usefulFiles)
    fileBase = usefulFiles{i};
    [w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,factors(i),phase_new);
end


%% Complex pulse
% folder = '../Daten/FROGS 10kHz_2mJ/4mJ/';
% folderContent = dir(folder);
% numberOfFiles = length(folderContent);
% usefulFiles = {};
% % Filtering of unwanted files necessary!
% fprintf('Searching for files with pattern: %s\n\n',pattern)
% for i=1:numberOfFiles
%     fileName = folderContent(i).name;
%     if ~isempty(regexp(fileName,pattern,'match'))
%         fileName = strrep(fileName,'.Speck.dat','');
%         fileName = strrep(fileName,'.Ek.dat','');
%         L = length(usefulFiles);
%         flag = true;
%         for j=1:L
%             if fileName==usefulFiles{j}
%                 flag = false;
%             end
%         end
%         if flag
%             usefulFiles{L+1} = fileName;
%             fprintf('File %d: %s\n',L+1,fileName)
%         end
%     end
% end
% if isempty(usefulFiles)
%     error('No files found that match the patter!')
% end
% for i=1:length(usefulFiles)
%     fileBase = usefulFiles{i};
%     [w_Sk,I_Sk,p_Sk,phase_new] = compensation_compensatePulse(folder,fileBase,1);
% end