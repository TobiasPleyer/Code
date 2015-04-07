clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
set(0,'DefaultFigureWindowStyle','docked') % Plots will be docked in the IDE main window
F_size_Label = 16;
colors = {'blue', 'red', 'green', 'cyan', 'magenta', 'black','yellow'};

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

figure()
[AX, ~, ~] = plotyy(ucmp_Speck_wavel,ucmp_Speck_int,ucmp_Speck_wavel,ucmp_Speck_phase);
title('The uncompressed pulse spectrum','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
figure()
[AX, ~, ~] = plotyy(cmp_Speck_wavel,cmp_Speck_int,cmp_Speck_wavel,cmp_Speck_phase);
title('The compressed pulse spectrum','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);

% In order to perform calculations we always need to bring all data on a
% common x-axis value basis
common_wavelength = (1012:0.1:1048)';
cmp_Speck_int = interp1(cmp_Speck_wavel,cmp_Speck_int,common_wavelength);
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
cmp_Speck_phase = interp1(cmp_Speck_wavel,cmp_Speck_phase,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);
theo_phase = interp1(theo_wavel,theo_phase,common_wavelength);

compressor_phase = cmp_Speck_phase - ucmp_Speck_phase;
phase_1bounce = compressor_phase ./ 20;
compressor_phase_one_more = cmp_Speck_phase + phase_1bounce; % What is the effect of one bounce more?
compressor_phase_one_less = cmp_Speck_phase - phase_1bounce; % What is the effect of one bounce less?
phase_after_theo = ucmp_Speck_phase - theo_phase;
phase_after_1M = ucmp_Speck_phase + phase_meas_1M;
phase_after_Full = ucmp_Speck_phase + phase_meas_Full;

%Calculate the compression as expected from the theoretical values and
%measurements
[Int_T,t_T,Ek_T] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,phase_after_theo,5);
[Int_1M,t_1M,Ek_1M] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,phase_after_1M,5);
[Int_Full,t_Full,Ek_Full] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,phase_after_Full,5);
[Int_1more,t_1more,Ek_1more] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,compressor_phase_one_more,5);
[Int_1less,t_1less,Ek_1less] = compressor_toTimeDomain(ucmp_Speck_int,common_wavelength,compressor_phase_one_less,5);

figure()
[AX, ~, ~] = plotyy(ucmp_Ek_time,ucmp_Ek_int,ucmp_Ek_time,ucmp_Ek_phase);
%title('The uncompressed pulse in time','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(AX(2),'YTick',-25:5:10)
set(AX(1),'Box','off')
set(AX,'FontSize',F_size_Label)

figure()
[AX, ~, ~] = plotyy(cmp_Ek_time,cmp_Ek_int,cmp_Ek_time,cmp_Ek_phase);
title('The compressed pulse in time after the compressor','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)

figure()
Int_meas = abs(trapz(cmp_Ek_time,cmp_Ek_int));
H1 = zeros(1,7);
H2 = zeros(1,6);
[AX, H1(1), H2(1)] = plotyy(cmp_Ek_time,cmp_Ek_int./Int_meas.*Int_F,cmp_Ek_time,cmp_Ek_phase);
set(H1(1),'Color',colors{1},'LineStyle','-')
set(H2(1),'Color',colors{1},'LineStyle','--')
hold(AX(1))
hold(AX(2))
[H1(2)] = plot(AX(1),t_F,abs(Ek_F).^2);
set(H1(2),'Color',colors{2},'LineStyle','-')
N = 210;
Ek_T = circshift(Ek_T',N);
[H1(3)] = plot(AX(1),t_T,abs(Ek_T).^2./Int_T.*Int_F);
[H2(2)] = plot(AX(2),t_T,-unwrap(angle(Ek_T)));
set(H1(3),'Color',colors{3},'LineStyle','-')
set(H2(2),'Color',colors{3},'LineStyle','--')
N = -22;
Ek_Full = circshift(Ek_Full',N);
[H1(4)] = plot(AX(1),t_Full,abs(Ek_Full).^2./Int_Full.*Int_F);
[H2(3)] = plot(AX(2),t_Full,-unwrap(angle(Ek_Full)));
set(H1(4),'Color',colors{4},'LineStyle','-')
set(H2(3),'Color',colors{4},'LineStyle','--')
N = -254;
Ek_1M = circshift(Ek_1M',N);
[H1(5)] = plot(AX(1),t_1M,abs(Ek_1M).^2./Int_1M.*Int_F);
[H2(4)] = plot(AX(2),t_1M,-unwrap(angle(Ek_1M)));
set(H1(5),'Color',colors{5},'LineStyle','-')
set(H2(4),'Color',colors{5},'LineStyle','--')
N = -12;
Ek_1more = circshift(Ek_1more',N);
[H1(6)] = plot(AX(1),t_1more,abs(Ek_1more).^2./Int_1more.*Int_F);
[H2(5)] = plot(AX(2),t_1more,-unwrap(angle(Ek_1more)));
set(H1(6),'Color',colors{6},'LineStyle','-')
set(H2(5),'Color',colors{6},'LineStyle','--')
N = 12;
Ek_1less = circshift(Ek_1less',N);
[H1(7)] = plot(AX(1),t_1less,abs(Ek_1less).^2./Int_1less.*Int_F);
[H2(6)] = plot(AX(2),t_1less,-unwrap(angle(Ek_1less)));
set(H1(7),'Color',[1,0.65,0.0],'LineStyle','-')
set(H2(6),'Color',[1,0.65,0.0],'LineStyle','--')
title(sprintf('Menlo FROG spectrum before compressor'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
xlim(AX(1),[-2000 2000])
xlim(AX(2),[-2000 2000])
% ylim(AX(2),[-100 150])
l = legend(H1,{'Measured', 'Fourier limit', 'Theoretical', 'Expected from FROG measurement', 'Expected from 1 Mirror', '1 bounce more', '1 bounce less'},'Location','NorthWest');

figure()
theo_phase = -theo_phase; % This is due to unfortunate sign conventions
theo_phase_tilted = compressor_alignCurve(common_wavelength,compressor_phase,theo_phase,50);
phase_meas_1M_tilted = compressor_alignCurve(common_wavelength,compressor_phase,phase_meas_1M,50);
phase_meas_Full_tilted = compressor_alignCurve(common_wavelength,compressor_phase,phase_meas_Full,50);
plot(common_wavelength, compressor_phase)
hold on
plot(common_wavelength, theo_phase_tilted,'r')
plot(common_wavelength, phase_meas_1M_tilted,'g')
plot(common_wavelength, phase_meas_Full_tilted,'c')
plot(common_wavelength, ucmp_Speck_int .* 35 - 25,'k')
hold off
xlim([1022 1038])
ylim([-25 10])
legend('Measured','Theoretical','FROG 1M','FROG Full')
title('Comparison of the calculated compressor phase, the theoretical phase and the phases determined beforehand','FontSize',F_size_Label-2)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)

figure()
phase_after_theo_tilted = compressor_alignCurve(common_wavelength,compressor_phase,phase_after_theo,50);
phase_after_1M_tilted = compressor_alignCurve(common_wavelength,compressor_phase,phase_after_1M,50);
phase_after_Full_tilted = compressor_alignCurve(common_wavelength,compressor_phase,phase_after_Full,50);
plot(common_wavelength, cmp_Speck_phase)
hold on
plot(common_wavelength, phase_after_theo_tilted,'r')
plot(common_wavelength, phase_after_1M_tilted,'g')
plot(common_wavelength, phase_after_Full_tilted,'c')
plot(common_wavelength, cmp_Speck_int .* 35 - 25,'k')
hold off
xlim([1022 1038])
ylim([-25 10])
legend('Measured','Theoretical','FROG 1M','FROG Full')
title(sprintf('Comparison of the measured phase after the compressor and the expected\nphases according to theory, and previous FROG measurements'),'FontSize',F_size_Label-2)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)