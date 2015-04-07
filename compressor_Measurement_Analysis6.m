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
cmp_Speck_int = interp1(cmp_Speck_wavel,cmp_Speck_int,common_wavelength);
ucmp_Speck_int = interp1(ucmp_Speck_wavel,ucmp_Speck_int,common_wavelength);
cmp_Speck_phase = interp1(cmp_Speck_wavel,cmp_Speck_phase,common_wavelength);
ucmp_Speck_phase = interp1(ucmp_Speck_wavel,ucmp_Speck_phase,common_wavelength);

compressor_phase = cmp_Speck_phase - ucmp_Speck_phase;
phase_1bounce = compressor_phase ./ 20;

bounces = [11 12 13];
n = 1;
N = length(bounces);
Specks = {};
for i=n:N
    Specks{i-n+1} = ucmp_Speck_phase + bounces(i).* phase_1bounce;
    [Int_N_less(i-n+1),t_N_less(:,i-n+1),Ek_N_less(:,i-n+1)] = compressor_toTimeDomain(cmp_Speck_int,common_wavelength,Specks{i-n+1},5);
end

Int_meas = abs(trapz(cmp_Ek_time,cmp_Ek_int));
figure()
plot(t_F,abs(Ek_F).^2,'Color',colors{1});
hold on
plot(t_C,abs(Ek_C).^2./Int_C.*Int_F,'Color',colors{2})
leg = {'Fourier limit',sprintf('Compression result last\ntime (20 bounces)')};
%Start to prepare the variables for the following plots
indices = {'003','004','006'};
expl = {'11 bounces @17.8 W','12 bounces @17.4 W','13 bounces @15.9 W'};
folder_base = 'T:/Tobias/Chirped Mirror Compressor/Analysis/20150326/Data/';
folder_cmp = folder_base;
for i=1:length(indices)
    filebase_cmp = [indices{i} '_FROG_Compressor_11bounces'];
    filename_cmp_Ek = [folder_cmp filebase_cmp Ek_suffix];
    filename_cmp_Speck = [folder_cmp filebase_cmp Speck_suffix];
    temp = dlmread(filename_cmp_Ek);
        Ek_time = temp(:,1);
        Ek_int = temp(:,2);
        Ek_phase = temp(:,3);
    temp = dlmread(filename_cmp_Speck);
        Speck_wavel = temp(:,1);
        Speck_int = temp(:,2);
        Speck_phase = temp(:,3);
    [Int,t,Ek] = compressor_toTimeDomain(Speck_int,Speck_wavel,Speck_phase,5);
    plot(t,abs(Ek).^2./Int.*Int_F,'Color',colors{2+i})
    leg{2+i} = expl{i};
end
hold off
h = legend(leg);
%title('Measured results for the possibility to improve the compressor','FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('intensity a.u.','FontSize',F_size_Label)
xlim([-2000 2000])
set(gca,'FontSize',F_size_Label)
set(h,'FontSize',10)
saveas(gcf,'../Bilder/Measured_Improvement','epsc')

%################### REFERENCE ############################################
% 
% figure()
% plot(t_F,abs(Ek_F).^2,'Color',colors{1});
% hold on
% plot(t_C,abs(Ek_C).^2./Int_C.*Int_F,'Color',colors{2})
% leg = {'Fourier limit','Compression result last time (20 bounces)'};
% %Start to repare the variables for the following plots
% indices = {'003','004','006'};
% expl = {'11 bounces @17.8 W','12 bounces @17.4 W','13 bounces @15.9 W'};
% folder_base = 'T:/Tobias/Chirped Mirror Compressor/Analysis/20150326/';
% folder_cmp = folder_base;
% for i=1:length(indices)
%     filebase_cmp = [indices{i} '_FROG_Compressor_11bounces'];
%     filename_cmp_Ek = [folder_cmp filebase_cmp Ek_suffix];
%     filename_cmp_Speck = [folder_cmp filebase_cmp Speck_suffix];
%     temp = dlmread(filename_cmp_Ek);
%         Ek_time = temp(:,1);
%         Ek_int = temp(:,2);
%         Ek_phase = temp(:,3);
%     temp = dlmread(filename_cmp_Speck);
%         Speck_wavel = temp(:,1);
%         Speck_int = temp(:,2);
%         Speck_phase = temp(:,3);
%     [Int,t,Ek] = compressor_toTimeDomain(Speck_int,Speck_wavel,Speck_phase,5);
%     [Int_F,~,~] = compressor_toTimeDomain(Speck_int,Speck_wavel,5);
%     plot(t,abs(Ek).^2./Int.*Int_F,'Color',colors{2+i})
%     leg{2+i} = expl{i};
% end
% hold off
% legend(leg)
% title('Measured results for the possibility to improve the compressor','FontSize',14)
% xlabel('time [fs]','FontSize',F_size_Label)
% ylabel('intensity a.u.','FontSize',F_size_Label)
% xlim([-2000 2000])
% set(gca,'FontSize',F_size_Label)

%################# END REFERENCE ##########################################

figure()
plot(t_F,abs(Ek_F).^2,'Color',colors{1});
hold on
%plot(cmp_Ek_time,cmp_Ek_int./Int_meas.*Int_F)
plot(t_C,abs(Ek_C).^2./Int_C.*Int_F,'Color',colors{2})
leg = {'Fourier limit',sprintf('Compression result last\ntime (20 bounces)')};
for i=n:N
    l = [num2str(bounces(i)) ' bounce(s)'];
    leg{2+i-n+1} = l;
    plot(t_N_less(:,i-n+1),abs(Ek_N_less(:,i-n+1)).^2./Int_N_less(i-n+1).*Int_F,'Color',colors{i-n+3});
end
h = legend(leg);
%title('Numerical calculations for the possibility to improve the compressor','FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('intensity a.u.','FontSize',F_size_Label)
xlim([-2000 2000])
set(gca,'FontSize',F_size_Label)
set(h,'FontSize',10)
saveas(gcf,'../Bilder/Calculated_Improvement','epsc')

HD536_files = {
    '003_FROG_Compressor_11bounces',
    '004_FROG_Compressor_11bounces',
    '006_FROG_Compressor_11bounces'
};
file_titles = {
    'Compressor_HD536_11_bounces_17,8W',
    'Compressor_HD536_12_bounces_17,4W',
    'Compressor_HD536_13_bounces_15,9W'
};

for i=1:length(HD536_files)
    auxiliary_plotPulse(HD536_files{i},true,file_titles{i});
end

%Observe if the spectrum changes significantly between 15.9 and 21.9 W
powers = zeros(1,11);
for i=1:11
    filename = sprintf('%03d_Power_Compressor_11bounces.txt',i);
    powers(i) = auxiliary_getPower(filename);
end
[m,i] = min(powers);
[M,I] = max(powers);
avg = mean(powers);
i_avg = auxiliary_find_closest_idx(powers,avg);
figure()
hold on
leg = {};
filename = sprintf('%03d_Spectrum_Compressor_11bounces.txt',i);
leg{1} = sprintf('Spectrum @%2.2f W',powers(i));
spectrum = dlmread(filename,'\t',[15,0,3600,1]);
plot(spectrum(:,1),spectrum(:,2)./max(spectrum(:,2)),'Color',colors{1})
filename = sprintf('%03d_Spectrum_Compressor_11bounces.txt',I);
leg{2} = sprintf('Spectrum @%2.2f W',powers(I));
spectrum = dlmread(filename,'\t',[15,0,3600,1]);
plot(spectrum(:,1),spectrum(:,2)./max(spectrum(:,2)),'Color',colors{2})
filename = sprintf('%03d_Spectrum_Compressor_11bounces.txt',i_avg);
leg{3} = sprintf('Spectrum @%2.2f W',powers(i_avg));
spectrum = dlmread(filename,'\t',[15,0,3600,1]);
plot(spectrum(:,1),spectrum(:,2)./max(spectrum(:,2)),'Color',colors{3})
%The 22W measurement
filename = '005_Spectrum_4bouncesHD535.txt';
leg{4} = 'Spectrum @21.88 W';
spectrum = dlmread(filename,'\t',[15,0,3600,1]);
plot(spectrum(:,1),spectrum(:,2)./max(spectrum(:,2)),'Color',colors{5})
hold off
h = legend(leg);
%title('Comparison between the spectra for different output powers','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('intensity a.u.','FontSize',F_size_Label)
xlim([1015 1045])
ylim([0 1.1])
set(gca,'FontSize',F_size_Label)
set(h,'FontSize',10)
saveas(gcf,'../Bilder/Comparison_Spectra','epsc')

names = {'006_FROG_Compressor_11bounces','007_FROG_Compressor_11bounces'};
expl = {'Compressor 13 HD536 @15.9 W','Compressor 12 HD536 + 1 HD535 @17.0 W'};
folder_base = 'T:/Tobias/Chirped Mirror Compressor/Analysis/20150326/Data/';
titl = 'Compression result when exchanging one HD536 mirror with a HD535 mirror';
auxiliary_plotPulses(folder_base,names,expl,titl)

HD535_files = {
    '002_FROG_8bouncesHD535',
    '003_FROG_4bouncesHD535',
    '005_FROG_4bouncesHD535'
};
file_titles = {
    'Compressor_HD535_8bounces_17,1W',
    'Compressor_HD535_4bounces_17,6W',
    'Compressor_HD535_4bounces_21,9W'
};

for i=1:length(HD535_files)
    auxiliary_plotPulse(HD535_files{i},true,file_titles{i});
end

%Find out the phase of the HD535 compressor
filename = '002_FROG_8bouncesHD535';
filename_Speck = [folder_base filename '.bin.Speck.dat'];
temp = dlmread(filename_Speck);
Speck_wavel = temp(:,1);
Speck_int = temp(:,2);
Speck_phase = temp(:,3);
Speck_int = interp1(Speck_wavel,Speck_int,common_wavelength);
Speck_phase = interp1(Speck_wavel,Speck_phase,common_wavelength);

phase_of_HD535 = Speck_phase - ucmp_Speck_phase;
figure()
[AX,~,~] = plotyy(common_wavelength,ucmp_Speck_int,common_wavelength,phase_of_HD535);
%title(sprintf('The phase of the compressor HD535 (8 bounces) drawn\nover the spectrum of the uncompressed pulse'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(AX,'FontSize',F_size_Label)
saveas(gcf,'../Bilder/Phase_HD535','epsc')

%Find out the phase of the HD536 compressor
filename = '004_FROG_Compressor_11bounces';
filename_Speck = [folder_base filename '.bin.Speck.dat'];
temp = dlmread(filename_Speck);
Speck_wavel = temp(:,1);
Speck_int = temp(:,2);
Speck_phase = temp(:,3);
Speck_int = interp1(Speck_wavel,Speck_int,common_wavelength);
Speck_phase = interp1(Speck_wavel,Speck_phase,common_wavelength);

phase_of_HD536 = Speck_phase - ucmp_Speck_phase;
figure()
[AX,~,~] = plotyy(common_wavelength,ucmp_Speck_int,common_wavelength,phase_of_HD536);
%title(sprintf('The phase of the compressor HD536 (12 bounces) drawn\nover the spectrum of the uncompressed pulse'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(AX,'FontSize',F_size_Label)
saveas(gcf,'../Bilder/Phase_HD536','epsc')