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
filename_narrow_sent = 'T:\Tobias\Chirped Mirror Compressor\Originals\007_FROG_AIR_Without_multipass_RTT3.645';

temp = dlmread(filename_ucmp_Ek);
ucmp_Ek_time = temp(:,1);
ucmp_Ek_int = temp(:,2);
ucmp_Ek_phase = temp(:,3);
temp = dlmread(filename_ucmp_Speck);
ucmp_Speck_wavel = temp(:,1);
ucmp_Speck_int = temp(:,2);
ucmp_Speck_phase = temp(:,3);
temp = dlmread([filename_narrow_sent '.bin.Ek.dat']);
sent_Ek_time = temp(:,1);
sent_Ek_int = temp(:,2);
sent_Ek_phase = temp(:,3);
temp = dlmread([filename_narrow_sent '.bin.Speck.dat']);
sent_Speck_wavel = temp(:,1);
sent_Speck_int = temp(:,2);
sent_Speck_phase = -temp(:,3);

clear temp;

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);

% In order to perform calculations we always need to bring all data on a
% common x-axis value basis
common_wavelength = (1015:0.1:1045)';
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
sent_Speck_int = interp1(sent_Speck_wavel,sent_Speck_int,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);
sent_Speck_phase = interp1(sent_Speck_wavel,sent_Speck_phase,common_wavelength);

[ucmp_P,ucmp_GD,ucmp_GDD] = compressor_GDD_from(common_wavelength,ucmp_Speck_phase,'+',1022,1038);
[sent_P,sent_GD,sent_GDD] = compressor_GDD_from(common_wavelength,sent_Speck_phase,'+',1022,1038);

% This is ugly code to remove the phase jump in the original
idx = auxiliary_find_closest_idx(common_wavelength,1029.9);
ucmp_Speck_phase2 = auxiliary_removeJump(ucmp_Speck_phase,idx,-1.9,-10);
    
phase_difference = ucmp_Speck_phase2 - sent_Speck_phase;

figure()
plot(common_wavelength,ucmp_Speck_phase)
hold on
plot(common_wavelength,ucmp_Speck_phase2,'g-')
plot(common_wavelength,sent_Speck_phase,'r-')
plot(common_wavelength,ucmp_Speck_int.*10-2,'b--')
plot(common_wavelength,sent_Speck_int.*10-2,'r--')
hold off
legend('Pulse sent into compressor',sprintf('Pulse sent into compressor\nwith removed phase jump'),'Pulse used for production')
title(sprintf('Direct comparison between the pulse used to design the\nmirrors and the pulse that was sent into the compressor'),'FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
xlim([1022 1038])
set(gca,'FontSize',F_size_Label)

fit_x = 1022:0.1:1038;
fit_y = interp1(common_wavelength,phase_difference,fit_x);
P = polyfit(fit_x,fit_y,5);
fit_y = polyval(P,fit_x);

figure()
[AX,~,~] = plotyy(common_wavelength,ucmp_Speck_int,common_wavelength,phase_difference);
hold(AX(2))
plot(AX(2),fit_x,fit_y,'r-')
title('The phase difference P_{sent in} - P_{production}','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
xlim(AX(1),[1022 1038])
xlim(AX(2),[1022 1038])
ylim(AX(2),[-5 0])
legend(AX(2),'phase difference','polynomial fit')
set(AX,'FontSize',F_size_Label)
set(AX(2),'YTick',-5:1:0)
set(AX(1),'Box','off')
set(AX,'FontSize',F_size_Label)

figure()
plot(common_wavelength,ucmp_GDD)
hold on
plot(common_wavelength,sent_GDD,'r-')
hold off
legend('Pulse sent into compressor','Pulse used for production')
title(sprintf('The GDDs of the measured uncompressed pulse and\nthe original pulse for whom the mirrors were designed'),'FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GDD [fs^2]','FontSize',F_size_Label)
xlim([1022 1038])
set(gca,'FontSize',F_size_Label)

[~,GDD] = compressor_HOD_from(fit_x,fit_y,2,10);
[~,TOD] = compressor_HOD_from(fit_x,fit_y,3,10);
[~,FOD] = compressor_HOD_from(fit_x,fit_y,4,10);

figure()
plot(fit_x,GDD,'r-')
title('The GDD for the phase difference P_{sent in} - P_{production}','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
set(gca,'FontSize',F_size_Label)
figure()
plot(fit_x,TOD,'r-')
title('The TOD for the phase difference P_{sent in} - P_{production}','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
set(gca,'FontSize',F_size_Label)
figure()
plot(fit_x,FOD,'r-')
title('The FOD for the phase difference P_{sent in} - P_{production}','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
set(gca,'FontSize',F_size_Label)