%##################### INIT #####################
clear

folder_base = 'T:\LEX_measurements\hybrid data\';
folder_path = '20150202\Tobi\';

cd([folder_base folder_path]) % change to the working directory
styles = {'r-', 'b-', 'g-', 'c-', 'm-', 'k-'};
colors = {'red', 'blue', 'green', 'cyan', 'magenta', 'black'};
F_size_Label = 22;
c = 299792458;

% Get the theoretical and measured values for the GD and GDD
Trubetskov = dlmread('Trubetskov_Mirror_Specification_Broad_Band.txt','\t',1,0);
Trubetskov_wavelength = Trubetskov(:,1);
Trubetskov_phase = Trubetskov(:,2);
Trubetskov_GD = Trubetskov(:,3);
Trubetskov_GDD = Trubetskov(:,4);

%##################### SPECTRUM #####################
%Load the spectrum data from file into variable
initial_spec = dlmread('Measured\001_SHG-FROG_Pure_Menlo_No_Mirrors_Reference_Spectrum.txt','\t',17,0);
initial_spec = initial_spec(:,1:2);
resultant_spec = dlmread('Measured\001_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo_HR4000_Reference_Spectrum.txt','\t',17,0);
resultant_spec = resultant_spec(:,1:2);

% Reduce the amount of unnecessary data (zeros)
lower_cutoff = 1000;
upper_cutoff = 1060;
initial_spec = initial_spec(initial_spec(:,1) > lower_cutoff & initial_spec(:,1) < upper_cutoff,:);
resultant_spec = resultant_spec(resultant_spec(:,1) > lower_cutoff & resultant_spec(:,1) < upper_cutoff,:);

%##################### FROG #####################
file_with_mirrors = '009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo';
file_without_mirrors = '004_SHG-FROG_Pure_Menlo_No_Mirrors';
file_end = '.bin.Speck.dat';
file_end2 = '.bin.Ek.dat';
file_folder = 'FROG/';

%##################### VALIDATION #####################
% To see that this is the correct procedure check the following call chain:
% T:\Software\FROG_v1.2.6_Lenard_backup\lib\frog\Frogger\FroggerSave.m
% T:\Software\FROG_v1.2.6_Lenard_backup\lib\io\esave.m
% T:\Software\FROG_v1.2.6_Lenard_backup\lib\Misc\IandP.m
% T:\Software\FROG_v1.2.6_Lenard_backup\lib\Misc\aandp.m
% Important is the line
%       Phase = -unwrap(angle(E));
% in the file aandp.m which comes from the sign convention of the phase:
%       S(w) = sqrt(I(w)) * exp(-iP(w))
% The minus does not come from the fft, so it is added and has to be
% accounted for in the ifft.
%--------------------------------------------------------------------------
% file_name = 'FROG/009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo.bin.Speck.dat';
% S = dlmread(file_name);
% wavelength_with_mirrors = S(:,1);
% spectral_intensity_with_mirrors = S(:,2);
% spectral_phase_with_mirrors = S(:,3);
% file_name = 'FROG/009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo.bin.Ek.dat';
% S = dlmread(file_name);
% time_with_mirrors = S(:,1);
% temporal_intensity_with_mirrors = S(:,2);
% temporal_phase_with_mirrors = S(:,3);
% figure()
% subplot(2,2,1)
% plotyy(wavelength_with_mirrors,spectral_intensity_with_mirrors,wavelength_with_mirrors,spectral_phase_with_mirrors)
% subplot(2,2,2)
% plotyy(time_with_mirrors,temporal_intensity_with_mirrors,time_with_mirrors,temporal_phase_with_mirrors)
% subplot(2,2,3)
% [Int,t,E] = compressor_toTimeDomain(spectral_intensity_with_mirrors,wavelength_with_mirrors,spectral_phase_with_mirrors);
% plotyy(t,abs(E).^2,t,unwrap(-angle(E)))
% subplot(2,2,4)
% plotyy(t,abs(E).^2-temporal_intensity_with_mirrors',t,unwrap(-angle(E))-temporal_phase_with_mirrors')


file_name = [file_folder file_with_mirrors file_end];
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
wavelength_with_mirrors = S(:,1);
spectral_intensity_with_mirrors = S(:,2);
spectral_phase_with_mirrors = -S(:,3); % The minus comes from the assumption that the mirrors are produced correctly, then the phase should have negative GDD
%spectral_phase_with_mirrors = -0.1*(wavelength_with_mirrors-1030).^2 + 0.01*(wavelength_with_mirrors-1030) - 5;

file_name = [file_folder file_without_mirrors file_end];
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
wavelength_without_mirrors = S(:,1);
spectral_intensity_without_mirrors = S(:,2);
spectral_phase_without_mirrors = -S(:,3); % The minus comes because we know the Menlo has positive chirp

%##################### MIRROR SIGN CHECK ####################
% Calculate the Fourier transforms of the pulse with measured mirror phase
% added and subtracted. Compare the results with the retrieved pulse to 
% understand if the produced mirrors implement a positive or negative chirp.
% We know that the Menlo laser is positivle chirped, i.e. the phase curve
% will have a positive second derivative in the spectral domain.
% Checking with the plots shows us that we have the correct phase for the
% Menlo retrieval.

figure()
subplot(2,2,1)
[AX, H1, H2] = plotyy(wavelength_without_mirrors,spectral_intensity_without_mirrors,wavelength_without_mirrors,spectral_phase_without_mirrors);
title(sprintf('Menlo FROG spectrum before compressor'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
xlim(AX(1),[1000 1060])
xlim(AX(2),[1000 1060])

subplot(2,2,2)
[AX, H1, H2] = plotyy(wavelength_with_mirrors,spectral_intensity_with_mirrors,wavelength_with_mirrors,spectral_phase_with_mirrors);
title(sprintf('Menlo FROG spectrum after compressor'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
xlim(AX(1),[1000 1060])
xlim(AX(2),[1000 1060])

subplot(2,2,3)
file_name = sprintf([file_folder file_with_mirrors file_end2],2);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
time_with_mirrors = S(:,1);
temporal_intensity_with_mirrors = fliplr(S(:,2)')';
temporal_phase_with_mirrors = -fliplr(S(:,3)')';
[AX, H1, H2] = plotyy(time_with_mirrors,temporal_intensity_with_mirrors,time_with_mirrors,temporal_phase_with_mirrors);
title(sprintf('Menlo FROG temporal measurement after compressor'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)

% Before we can start to find the correct phase difference we have to make
% sure that the indices offset in both files is identical. This is not the
% case:
% >> auxiliary_find_closest_idx(wavelength_with_mirrors,1049)
%       ans = 167
% >> auxiliary_find_closest_idx(wavelength_without_mirrors,1049)
%       ans = 204
% To assure that we interpolate the data on a common X-scale

common_wavelength = 1000:0.1:1060;
spectral_intensity_without_mirrors = interp1(wavelength_without_mirrors,spectral_intensity_without_mirrors,common_wavelength);
spectral_intensity_with_mirrors = interp1(wavelength_with_mirrors,spectral_intensity_with_mirrors,common_wavelength);
spectral_phase_without_mirrors = interp1(wavelength_without_mirrors,spectral_phase_without_mirrors,common_wavelength);
spectral_phase_with_mirrors = interp1(wavelength_with_mirrors,spectral_phase_with_mirrors,common_wavelength);
spectral_phase_diff = spectral_phase_with_mirrors-spectral_phase_without_mirrors;
spectral_phase_diff2 = -spectral_phase_with_mirrors-spectral_phase_without_mirrors;

subplot(2,2,4)
% We don't know if we measured the real compressor result pulse or the time
% flipped version of it. So we end up with four combinations we have to
% check in order to find out what the phase of our chirped mirror is.
%--------------------------------------------------------------------------
Phi_mirror1 = spectral_phase_without_mirrors + spectral_phase_diff;
Phi_mirror2 = spectral_phase_without_mirrors - spectral_phase_diff;
Phi_mirror3 = spectral_phase_without_mirrors + spectral_phase_diff2;
Phi_mirror4 = spectral_phase_without_mirrors - spectral_phase_diff2;

[Int_1,t_1,Ek_1] = compressor_toTimeDomain(spectral_intensity_without_mirrors,common_wavelength,Phi_mirror1);
[Int_2,t_2,Ek_2] = compressor_toTimeDomain(spectral_intensity_without_mirrors,common_wavelength,Phi_mirror2);
[Int_3,t_3,Ek_3] = compressor_toTimeDomain(spectral_intensity_without_mirrors,common_wavelength,Phi_mirror3);
[Int_4,t_4,Ek_4] = compressor_toTimeDomain(spectral_intensity_without_mirrors,common_wavelength,Phi_mirror4);
Int_FROG = abs(trapz(time_with_mirrors,temporal_intensity_with_mirrors));

H1 = zeros(1,5); % Allocate the array for the handles
H2 = zeros(1,5);
[AX,H1(1),H2(1)] = plotyy(t_1,abs(Ek_1).^2./Int_1,t_1,unwrap(-angle(Ek_1)));
title('Fourier transformations of several phase combinations','FontSize',14)
set(AX(1),'YColor','black')
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','Intensity/Integral','FontSize',F_size_Label)
set(AX(2),'YColor','black')
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(H1(1),'Color',colors{1},'LineWidth',2)
set(H2(1),'Color',colors{1},'LineStyle','--')
hold(AX(1))
hold(AX(2))
[H1(2)] = plot(AX(1),t_2,abs(Ek_2).^2./Int_2);
[H2(2)] = plot(AX(2),t_2,unwrap(-angle(Ek_2)));
set(H1(2),'Color',colors{2},'LineWidth',2)
set(H2(2),'Color',colors{2},'LineStyle','--')
[H1(3)] = plot(AX(1),t_3,abs(Ek_3).^2./Int_3);
[H2(3)] = plot(AX(2),t_3,unwrap(-angle(Ek_3)));
set(H1(3),'Color',colors{3},'LineWidth',2)
set(H2(3),'Color',colors{3},'LineStyle','--')
[H1(4)] = plot(AX(1),t_4,abs(Ek_4).^2./Int_4);
[H2(4)] = plot(AX(2),t_4,unwrap(-angle(Ek_4)));
set(H1(4),'Color',colors{4},'LineWidth',2)
set(H2(4),'Color',colors{4},'LineStyle','--')
[H1(5)] = plot(AX(1),time_with_mirrors,temporal_intensity_with_mirrors./Int_FROG);
[H2(5)] = plot(AX(2),time_with_mirrors,temporal_phase_with_mirrors);
set(H1(5),'Color',colors{5},'LineWidth',2)
set(H2(5),'Color',colors{5},'LineStyle','--')
ylim(AX(1),[0 0.0013])
xlim(AX(1),[-10000 10000])
xlim(AX(2),[-10000 10000])
l = legend(H1,{'V1', 'V2', 'V3', 'V4', 'FROG measurement'},'Location','NorthWest');
set(l,'FontSize',F_size_Label);


%##################### GD/GDD ####################
poly = polyfit(common_wavelength,spectral_phase_diff,6);
phase_diff_poly_fit = polyval(poly,common_wavelength);
% We know the phase in wavelength domain -> transform to frequency space
% and calculate the GD and GDD
common_frequency = (2*pi*c) ./ (common_wavelength .* 1e-9);
poly = polyfit(common_frequency,spectral_phase_diff,6);
poly_deriv = polyder(poly);
poly_GD = polyval(poly_deriv,common_frequency) .* 1e15;
poly_deriv = polyder(poly_deriv);
poly_GDD = polyval(poly_deriv,common_frequency) .* 1e30;

figure()
subplot(2,2,1)
leg = {};
plot(common_wavelength,spectral_phase_diff,'r-','LineWidth',2)
leg{1} = sprintf('Measured Phase');
hold on
plot(common_wavelength,phase_diff_poly_fit,'g-','LineWidth',2)
leg{2} = sprintf('Polynomial Fit');
hold off
xlim([1000 1060])
title(sprintf('Measured and fitted Phase Function of the\nBroad Band Chirped Mirror Compressor'),'FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
l = legend(leg,'Location','NorthWest');
set(l,'FontSize',F_size_Label);

subplot(2,2,2)
% Before we can plot the data we have to make sure that they have the same
% offset from the axis in order to be comparable. The value of the
% theoretical measurement is taken as zero point
idx_end_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1050);
idx_end_M = auxiliary_find_closest_idx(common_wavelength,1050);
d_M = Trubetskov_GD(idx_end_T) - poly_GD(idx_end_M);
plot(Trubetskov_wavelength,Trubetskov_GD,'r-')
hold on
plot(common_wavelength,poly_GD+d_M,'g-')
hold off
xlim([1000 1060])
title(sprintf('Comparison of theoretical and measured GD\nvalues for the chirped mirror compressor'),'FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GD [fs]','FontSize',F_size_Label)
l = legend('Calculated by M. Trubetskov','Measured by T. Pleyer','Location','SouthEast');
set(l,'FontSize',F_size_Label);
subplot(2,2,3)
plot(Trubetskov_wavelength,Trubetskov_GDD,'r-')
hold on
plot(common_wavelength,poly_GDD,'g-')
hold off
xlim([1000 1060])
title(sprintf('Comparison of theoretical and measured\nGDD values for the chirped mirrors'),'FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GDD [fs^2]','FontSize',F_size_Label)
l = legend('Calculated by M. Trubetskov','Measured by T. Pleyer','Location','NorthWest');
set(l,'FontSize',F_size_Label);