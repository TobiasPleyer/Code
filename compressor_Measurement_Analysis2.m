clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
set(0,'DefaultFigureWindowStyle','docked') % Plots will be docked in the IDE main window
F_size_Label = 16;
colors = {'blue', 'red', 'green', 'cyan', 'magenta', 'black', 'yellow'};
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

%Calculate the Fourier limit for the uncompressed pulse
[Int_F,t_F,Ek_F] = compressor_toTimeDomain(ucmp_Speck_int,ucmp_Speck_wavel,5);
[Int_C,t_C,Ek_C] = compressor_toTimeDomain(cmp_Speck_int,cmp_Speck_wavel,cmp_Speck_phase,5);

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

bounces = [8 16 20];
n = 1;
N = length(bounces);
Specks = {};
for i=n:N
    Specks{i-n+1} = cmp_Speck_phase - bounces(i).* phase_1bounce;
    [Int_N_less(i-n+1),t_N_less(:,i-n+1),Ek_N_less(:,i-n+1)] = compressor_toTimeDomain(cmp_Speck_int,common_wavelength,Specks{i-n+1},5);
end

figure()
[AX, ~, ~] = plotyy(ucmp_Ek_time,ucmp_Ek_int,ucmp_Ek_time,ucmp_Ek_phase);
[t_min,t_max,Int] = auxiliary_findFWHM(ucmp_Ek_time,ucmp_Ek_int);
%title('The uncompressed pulse in time','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(AX,'FontSize',16)
set(AX(1),'Box','off')
ylim(AX(2),[-20 2])
set(AX(2),'YTick',-20:5:2)
xlim(AX(1),[-3000 3000])
xlim(AX(2),[-3000 3000])
hold(AX(1))
plot(AX(1),[t_min,t_min,t_max,t_max],[0,0.5,0.5,0],'r-')
label = sprintf('FWHM: %2.0f fs',t_max-t_min);
text((t_max+t_min)/2-length(label)*50,0.45,label,'FontSize',16)
saveas(gcf,'../Bilder/Additional/Uncompressed','epsc')

figure()
[Int,t,Ek] = compressor_toTimeDomain(cmp_Speck_int,common_wavelength,cmp_Speck_phase,5);
[t_min,t_max,Int] = auxiliary_findFWHM(t,abs(Ek).^2);
%[AX, ~, ~] = plotyy(cmp_Ek_time,cmp_Ek_int,cmp_Ek_time,cmp_Ek_phase);
[AX, ~, ~] = plotyy(t,abs(Ek).^2,cmp_Ek_time,cmp_Ek_phase);
%title('The compressed pulse in time after the compressor','FontSize',F_size_Label)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(AX,'FontSize',16)
set(AX(1),'Box','off')
ylim(AX(2),[-5 30])
set(AX(2),'YTick',-5:5:30)
xlim(AX(1),[-3000 3000])
xlim(AX(2),[-3000 3000])
hold(AX(1))
plot(AX(1),[t_min,t_min,t_max,t_max],[0,0.5,0.5,0],'r-')
label = sprintf('FWHM: %2.0f fs',t_max-t_min);
text(t_max + 40,0.45,label,'FontSize',16)
saveas(gcf,'../Bilder/Additional/Compressed_20Bounces','epsc')

figure()
plot(t_C,abs(Ek_C).^2./Int_C.*Int_F)
hold on
plot(t_F,abs(Ek_F).^2,'Color',colors{2});
leg  = {'Measured','Fourier limit'};
for i=n:N
    l = [num2str(bounces(i)) ' bounce(s) less'];
    leg{2+i-n+1} = l;
    plot(t_N_less(:,i-n+1),abs(Ek_N_less(:,i-n+1)).^2./Int_N_less(i-n+1).*Int_F,'Color',colors{mod(2+i-n,8)+1});
end
%title(sprintf('Compression result for reduced number of bounces'),'FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('intensity [a.u.]','FontSize',F_size_Label)
xlim([-2000 2000])
ylim([0 1.1])
l = legend(leg,'Location','NorthWest');
set(gca,'FontSize',F_size_Label)
set(l,'FontSize',10)
saveas(gcf,'../Bilder/Additional/Optimization_20Bounces','epsc')

figure()
plot(common_wavelength,cmp_Speck_phase,'Color',colors{1})
leg = {'Compressed with current setup'};
hold on
for i=n:N
    if bounces(i)==8
        l = [num2str(bounces(i)) ' bounce(s) less (optimum)'];
    elseif bounces(i)==20
        l = [num2str(bounces(i)) ' bounce(s) less (uncompressed case)'];
    else
        l = [num2str(bounces(i)) ' bounce(s) less'];
    end
    leg{i-n+2} = l;
    plot(common_wavelength,Specks{i-n+1},'Color',colors{mod(1+i-n,8)+1})
end
plot(common_wavelength,cmp_Speck_int.*16-10,'k-')
hold off
xlim([1022 1038])
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
%title('View of the spectra for configurations with different numbers of mirros','FontSize',14)
l = legend(leg,'Location','NorthEast');
set(gca,'FontSize',F_size_Label)
set(l,'FontSize',10)
saveas(gcf,'../Bilder/Additional/Optimization_20Bounces_Difference','epsc')

[t_min,t_max,Int] = auxiliary_findFWHM(t_C,abs(Ek_C).^2);
fprintf('FWHM with current setup: %f\n',t_max-t_min)
[t_min,t_max,Int] = auxiliary_findFWHM(t_N_less(:,1),abs(Ek_N_less(:,1)).^2);
fprintf('FWHM for optimum setup: %f\n',t_max-t_min)
[t_min,t_max,Int] = auxiliary_findFWHM(t_F,abs(Ek_F).^2);
fprintf('FWHM of Fourier limit: %f\n',t_max-t_min)