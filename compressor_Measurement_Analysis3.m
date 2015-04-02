clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
set(0,'DefaultFigureWindowStyle','docked') % Plots will be docked in the IDE main window
F_size_Label = 16;

folder_base = 'T:/Tobias/Chirped Mirror Compressor/Originals/Frogs/';
folder_base2 = 'T:\Tobias\Chirped Mirror Compressor\Analysis\20150326\Data\';
index = '001';
index2 = '003';
folder_cmp = [folder_base 'Compressed/' index '/'];
folder_cmp2 = folder_base2;
folder_ucmp = [folder_base 'Uncompressed/' index '/'];

filebase_cmp = [index '_FROG_Compressor1'];
filebase_cmp2 = [index2 '_FROG_Compressor_11bounces'];
filebase_ucmp = [index '_FROG_Test2'];

Ek_suffix = '.bin.Ek.dat';
Speck_suffix = '.bin.Speck.dat';

filename_cmp_Ek = [folder_cmp filebase_cmp Ek_suffix];
filename_cmp_Ek2 = [folder_cmp2 filebase_cmp2 Ek_suffix];
filename_ucmp_Ek = [folder_ucmp filebase_ucmp Ek_suffix];
filename_cmp_Speck = [folder_cmp filebase_cmp Speck_suffix];
filename_cmp_Speck2 = [folder_cmp2 filebase_cmp2 Speck_suffix];
filename_ucmp_Speck = [folder_ucmp filebase_ucmp Speck_suffix];
filename_narrow_Trubetskov  = 'T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Narrow_Band.txt';

% The phase data for the 1 mirror case and full compressor is saved in this
% .mat file. Be aware that the 1 mirror phase is already multiplied by the
% number of mirrors in the compressor.
load('T:/Tobias/Chirped Mirror Compressor/Analysis/Narrow_Variables_Comparison_1M_Full')

% I am still not sure which is the correct spectral phase difference that
% results from the measurements. For now I will choose the _v2 variables.
phase_meas_Full = narrow_Full_v2;
phase_meas_1M = narrow_1M_v2;

temp = dlmread(filename_cmp_Ek);
cmp_Ek_time = temp(:,1);
cmp_Ek_int = temp(:,2);
cmp_Ek_phase = temp(:,3);
temp = dlmread(filename_cmp_Speck);
cmp_Speck_wavel = temp(:,1);
cmp_Speck_int = temp(:,2);
cmp_Speck_phase = temp(:,3);
temp = dlmread(filename_cmp_Ek2);
cmp_Ek_time2 = temp(:,1);
cmp_Ek_int2 = temp(:,2);
cmp_Ek_phase2 = temp(:,3);
temp = dlmread(filename_cmp_Speck2);
cmp_Speck_wavel2 = temp(:,1);
cmp_Speck_int2 = temp(:,2);
cmp_Speck_phase2 = temp(:,3);
temp = dlmread(filename_ucmp_Ek);
ucmp_Ek_time = temp(:,1);
ucmp_Ek_int = temp(:,2);
ucmp_Ek_phase = temp(:,3);
temp = dlmread(filename_ucmp_Speck);
ucmp_Speck_wavel = temp(:,1);
ucmp_Speck_int = temp(:,2);
ucmp_Speck_phase = temp(:,3);
temp = dlmread(filename_narrow_Trubetskov,'\t',1,0);
theo_wavel = temp(:,1);
theo_phase = 20.*temp(:,2);

clear temp;

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);
[Int_C,t_C,Ek_C] = compressor_toTimeDomain(cmp_Speck_int,cmp_Speck_wavel,cmp_Speck_phase,5);
[Int_C2,t_C2,Ek_C2] = compressor_toTimeDomain(cmp_Speck_int2,cmp_Speck_wavel2,cmp_Speck_phase2,5);

% In order to perform calculations we always need to bring all data on a
% common x-axis value basis
common_wavelength = (1012:0.1:1048)';
cmp_Speck_int = interp1(cmp_Speck_wavel,cmp_Speck_int,common_wavelength);
cmp_Speck_int2 = interp1(cmp_Speck_wavel2,cmp_Speck_int2,common_wavelength);
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
cmp_Speck_phase = interp1(cmp_Speck_wavel,cmp_Speck_phase,common_wavelength);
cmp_Speck_phase2 = interp1(cmp_Speck_wavel2,cmp_Speck_phase2,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);
theo_phase = interp1(theo_wavel,theo_phase,common_wavelength);

compressor_phase = cmp_Speck_phase - ucmp_Speck_phase;
phase_1bounce = compressor_phase ./ 20;

[P,GD,GDD] = compressor_GDD_from(common_wavelength,cmp_Speck_phase);
figure()
plot(common_wavelength,GDD)
title('The GDD of the compressed pulse')
xlabel('wavelength [nm]')
ylabel('GDD [fs^2]')
xlim([1015 1045])

%See if the pre-pulse comes from the fact that we drop off to the sides
w = 1026;
W = 1034;
idx1 = auxiliary_find_closest_idx(common_wavelength,w);
idx2 = auxiliary_find_closest_idx(common_wavelength,W);
val1 = cmp_Speck_phase(idx1);
val2 = cmp_Speck_phase(idx2);
Speck_no_drop = cmp_Speck_phase;
Speck_no_drop(1:idx1) = val1;
Speck_no_drop(idx2:end) = val2;

[Int_no_drop,t_no_drop,Ek_no_drop] = compressor_toTimeDomain(cmp_Speck_int,common_wavelength,Speck_no_drop,5);

plot(common_wavelength,cmp_Speck_phase,'LineWidth',2)
hold on
plot(common_wavelength,Speck_no_drop,'r-')
plot(common_wavelength,cmp_Speck_int.*10-9,'k-')
hold off
legend('Measured phase','Modified phase')
title('Preventing the phase from dropping off at the edges','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
xlim([1022 1038])
set(gca,'FontSize',F_size_Label)

Int_meas = abs(trapz(cmp_Ek_time,cmp_Ek_int));
figure()
plot(t_F,abs(Ek_F).^2,'g-');
hold on
%plot(cmp_Ek_time,cmp_Ek_int./Int_meas.*Int_F)
plot(t_C,abs(Ek_C).^2./Int_C.*Int_F)
plot(t_C2,abs(Ek_C2).^2./Int_C2.*Int_F,'m-')
%plot(t_no_drop,abs(Ek_no_drop).^2./Int_no_drop.*Int_F,'r-')
hold off
legend('Fourier limit','20 bounces','11 bounces')
title('Preventing the phase from dropping off at the edges (resulting from over compensation)','FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('intensity a.u.','FontSize',F_size_Label)
xlim([-2000 2000])
set(gca,'FontSize',F_size_Label)