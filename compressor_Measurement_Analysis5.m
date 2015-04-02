clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
set(0,'DefaultFigureWindowStyle','docked') % Plots will be docked in the IDE main window
F_size_Label = 16;

folder_base = 'T:/Tobias/Chirped Mirror Compressor/Originals/Frogs/';
index = '001';
folder_cmp = [folder_base 'Compressed/' index '/'];
folder_ucmp = [folder_base 'Uncompressed/' index '/'];

filebase_cmp = [index '_FROG_Compressor1'];
filebase_ucmp = [index '_FROG_Test2'];

Ek_suffix = '.bin.Ek.dat';
Speck_suffix = '.bin.Speck.dat';

filename_ucmp_Ek = [folder_ucmp filebase_ucmp Ek_suffix];
filename_ucmp_Speck = [folder_ucmp filebase_ucmp Speck_suffix];

temp = dlmread(filename_ucmp_Ek);
ucmp_Ek_time = temp(:,1);
ucmp_Ek_int = temp(:,2);
ucmp_Ek_phase = temp(:,3);
temp = dlmread(filename_ucmp_Speck);
ucmp_Speck_wavel = temp(:,1);
ucmp_Speck_int = temp(:,2);
ucmp_Speck_phase = temp(:,3);

% The phase data for the 1 mirror case and full compressor is saved in this
% .mat file. Be aware that the 1 mirror phase is already multiplied by the
% number of mirrors in the compressor.
load('T:/Tobias/Chirped Mirror Compressor/Analysis/Broad_Variables_Comparison_1M_Full')

%Bad naming convention...
phase_meas_Full = P_full;
phase_meas_1M = P_1M;

clear temp;

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);

% In order to perform calculations we always need to bring all data on a
% common x-axis value basis
common_wavelength = (1015:0.1:1045)';
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);
P_full = interp1(wavelength,P_full,common_wavelength);

broad_compressed_phase = ucmp_Speck_phase-P_full;
[~,t_broad,Ek_broad] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,broad_compressed_phase);
[Int_broad,t_broad2,Ek_broad2] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,broad_compressed_phase,5);

figure()
P_full = compressor_alignCurve(common_wavelength,ucmp_Speck_phase,P_full);
H2 = zeros(1,2);
[AX,H1,H2(1)] = plotyy(common_wavelength,ucmp_Speck_int,common_wavelength,ucmp_Speck_phase);
hold(AX(2))
[H(2)] = plot(AX(2),common_wavelength,P_full,'r-');

figure()
[AX,~,~] = plotyy(common_wavelength,ucmp_Speck_int,common_wavelength,broad_compressed_phase);

[t_min,t_max,Int] = auxiliary_findFWHM(t_broad2,abs(Ek_broad2).^2);
figure()
[AX,~,~] = plotyy(t_broad2,abs(Ek_broad2).^2,t_broad,-unwrap(angle(Ek_broad)));
hold(AX(1))
plot(AX(1),[t_min,t_min,t_max,t_max],[0,0.5,0.5,0],'r-')
xlim(AX(1),[-3000 2000])
xlim(AX(2),[-3000 2000])

