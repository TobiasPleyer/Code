clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
set(0,'DefaultFigureWindowStyle','docked') % Plots will be docked in the IDE main window
F_size_Label = 16;
colors = {[0 0 1], [1 0 0], [0 1 0], [0 1 1], [1 0 1], [0 0 0], [1 0.6 0]};

folder_base = 'T:/Tobias/Chirped Mirror Compressor/Originals/Frogs/';
index = '001';
folder_cmp = [folder_base 'Compressed/' index '/'];
folder_ucmp = [folder_base 'Uncompressed/' index '/'];

filebase_cmp = [index '_FROG_Compressor1'];
filebase_ucmp = [index '_FROG_Test2'];

Ek_suffix = '.bin.Ek.dat';
Speck_suffix = '.bin.Speck.dat';

filename_cmp_Ek = [folder_cmp filebase_cmp Ek_suffix];
filename_ucmp_Ek = [folder_ucmp filebase_ucmp Ek_suffix];
filename_cmp_Speck = [folder_cmp filebase_cmp Speck_suffix];
filename_ucmp_Speck = [folder_ucmp filebase_ucmp Speck_suffix];

temp = dlmread(filename_cmp_Ek);
cmp_Ek_time = temp(:,1);
cmp_Ek_int = temp(:,2);
cmp_Ek_phase = temp(:,3);
temp = dlmread(filename_cmp_Speck);
cmp_Speck_wavel = temp(:,1);
cmp_Speck_int = temp(:,2);
cmp_Speck_phase = temp(:,3);
temp = dlmread(filename_ucmp_Ek);
ucmp_Ek_time = temp(:,1);
ucmp_Ek_int = temp(:,2);
ucmp_Ek_phase = temp(:,3);
temp = dlmread(filename_ucmp_Speck);
ucmp_Speck_wavel = temp(:,1);
ucmp_Speck_int = temp(:,2);
ucmp_Speck_phase = temp(:,3);

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);
[Int_C,t_C,Ek_C] = compressor_toTimeDomain(cmp_Speck_int,cmp_Speck_wavel,cmp_Speck_phase,5);

% In order to perform calculations we always need to bring all data on a
% common x-axis value basis
common_wavelength = (1012:0.1:1048)';
common_wavelength2 = (1015:0.1:1045)';
cmp_Speck_int = interp1(cmp_Speck_wavel,cmp_Speck_int,common_wavelength);
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
cmp_Speck_phase = interp1(cmp_Speck_wavel,cmp_Speck_phase,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);

compressor_phase = cmp_Speck_phase - ucmp_Speck_phase;
phase_1bounce = compressor_phase ./ 20;

load('T:\Tobias\Chirped Mirror Compressor\Analysis\Sent_Ucmp_Diff.mat')

HD536_files = {
    '003_FROG_Compressor_11bounces',
    '004_FROG_Compressor_11bounces',
    '006_FROG_Compressor_11bounces'
};

for i=1:3
    filename = [HD536_files{i} '.bin.Speck.dat'];
    temp = dlmread(filename);
    Speck_wavel = temp(:,1);
    Speck_int = temp(:,2);
    Speck_phase = temp(:,3);
    Speck_int = interp1(Speck_wavel,Speck_int,common_wavelength);
    Speck_phase = interp1(Speck_wavel,Speck_phase,common_wavelength);
    %We want to tilt Sent_Ucmp_Diff so that it is better comparable with the
    %residual phase of the pulse
    phase = interp1(common_wavelength,Speck_phase,common_wavelength2);
    Sent_Ucmp_Diff = auxiliary_tiltCurve(common_wavelength2,phase,Sent_Ucmp_Diff,1030,80);
    figure()
    [AX,H1,H2] = plotyy(common_wavelength,Speck_int,common_wavelength,Speck_phase);
    hold(AX(2))
    plot(AX(2),(1015:0.1:1045)',Sent_Ucmp_Diff,'r-')
    set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',16)
    set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',16)
    set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',16)
    set(AX(1),'Box','off')
    xlim(AX(1),[1020 1040])
    xlim(AX(2),[1020 1040])
    ylim(AX(2),[-5 10])
    set(AX(2),'YTick',-5:4:10)
    l = legend(AX(2),sprintf('Compressed phase\nafter %d bounces',10+i),'P_s - P_d','Location','NorthEast');
    set(AX,'FontSize',16)
    set(l,'FontSize',12)
    saveas(gcf,sprintf('T:/Tobias/Chirped Mirror Compressor/Analysis/20150326/Bilder/Additional/Comp_Res_Diff%d_bounces',10+i),'epsc')
end