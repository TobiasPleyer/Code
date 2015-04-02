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
% Vova = dlmread('Vova_Mirror_Data_Narrow_Band.dat','\t',2,0);
% Vova_wavelength = Vova(:,1);
% Vova_GD = Vova(:,2);
% Vova_GDD = Vova(:,3);


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

figure()

%##################### FROG #####################
with_mirrors = '_SHG-FROG_Compressor_angle_check_Reference_M1';
without_mirrors = '_SHG-FROG_Menlo_no_mirrors';
file_counter = '%03d';
file_end = '.bin.Speck.dat';
file_folder = 'FROG/';



%##################### subplot(2,2,1) #####################
N = 2; % Number of measurements made
fftshift_indices = [1]; % These indicate which plots need an fftshift applied
leg = {}; % the legend object to be displayed
% file_name = sprintf([file_folder file_counter with_mirrors file_end],1);

disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
leg{1} = sprintf(file_counter,1);
wavelength = S(:,1);
intensity = S(:,2);
phase = S(:,3);
H1 = zeros(1,N); % Allocate the array for the handles
H2 = zeros(1,N);
wavelengths_with_mirrors = zeros(256,N);
intensities_with_mirrors = zeros(256,N);
phases_with_mirrors = zeros(256,N);
if any(fftshift_indices == 1) > 0
    intensities_with_mirrors(:,1) = fftshift(intensity);
    wavelengths_with_mirrors(:,1) = wavelength;
    phases_with_mirrors(:,1) = unwrap(fftshift(phase));
else
    intensities_with_mirrors(:,1) = intensity;
    wavelengths_with_mirrors(:,1) = wavelength;
    phases_with_mirrors(:,1) = phase;
end
[AX,H1(1),H2(1)] = plotyy(wavelengths_with_mirrors(:,1),intensities_with_mirrors(:,1),wavelengths_with_mirrors(:,1),phases_with_mirrors(:,1));
title(sprintf('FROG retrieval of the pulse with\ntwo bounces off the chirped mirror'),'FontSize',14)
set(AX(1),'YColor','black')
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(AX(2),'YColor','black')
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(H1(1),'Color',colors{1},'LineWidth',2)
set(H2(1),'Color',colors{1},'LineStyle','--')
hold(AX(1))
hold(AX(2))
for n=2:N
    file_name = sprintf([file_folder file_counter with_mirrors file_end],n);
    disp(file_name)
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    leg{n} = sprintf(file_counter,n);
    wavelength = S(:,1);
    intensity = S(:,2);
    phase = S(:,3);
    if any(fftshift_indices == n) > 0
        intensities_with_mirrors(:,n) = fftshift(intensity);
        wavelengths_with_mirrors(:,n) = wavelength;
        phases_with_mirrors(:,n) = unwrap(fftshift(phase));
    else
        intensities_with_mirrors(:,n) = intensity;
        wavelengths_with_mirrors(:,n) = wavelength;
        phases_with_mirrors(:,n) = -1.*phase;
    end
    [H1(n)] = plot(AX(1),wavelengths_with_mirrors(:,n),intensities_with_mirrors(:,n));
    [H2(n)] = plot(AX(2),wavelengths_with_mirrors(:,n),phases_with_mirrors(:,n));
    set(H1(n),'Color',colors{n},'LineWidth',2)
    set(H2(n),'Color',colors{n},'LineStyle','--')
end
xlim(AX(1),[1000 1060])
xlim(AX(2),[1000 1060])
l = legend(H1,leg,'Location','South');
set(l,'FontSize',F_size_Label);
t = get(l,'Title');
set(t,'String','Measurement')
%legend([H1 H2],['Test1';'Test2'],['Muh';'Mäh'])



%##################### subplot(2,2,2) ####################
N = 2; % Number of measurements made
figure()
fftshift_indices = [1 2]; % These indicate which plots need an fftshift applied
%leg = {}; % the legend object to be displayed
file_name = sprintf([file_folder file_counter without_mirrors file_end],1);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
leg{1} = sprintf(file_counter,1);
wavelength = S(:,1);
intensity = S(:,2);
phase = S(:,3);
H1 = zeros(1,N); % Allocate the array for the handles
H2 = zeros(1,N);
wavelengths_without_mirrors = zeros(256,N);
intensities_without_mirrors = zeros(256,N);
phases_without_mirrors = zeros(256,N);
if any(fftshift_indices == 1) > 0
    intensities_without_mirrors(:,1) = fftshift(intensity);
    wavelengths_without_mirrors(:,1) = wavelength;
    phases_without_mirrors(:,1) = unwrap(-1.*fftshift(phase));
else
    intensities_without_mirrors(:,1) = intensity;
    wavelengths_without_mirrors(:,1) = wavelength;
    phases_without_mirrors(:,1) = phase;
end
[AX,H1(1),H2(1)] = plotyy(wavelengths_without_mirrors(:,1),intensities_without_mirrors(:,1),wavelengths_without_mirrors(:,1),phases_without_mirrors(:,1));
title(sprintf('FROG retrieval of the pulse directly out of the Menlo'),'FontSize',14)
set(AX(1),'YColor','black')
set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(AX(2),'YColor','black')
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
set(H1(1),'Color',colors{1},'LineWidth',2)
set(H2(1),'Color',colors{1},'LineStyle','--')
hold(AX(1))
hold(AX(2))
for n=2:N
    file_name = sprintf([file_folder file_counter without_mirrors file_end],n);
    disp(file_name)
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    leg{n} = sprintf(file_counter,n);
    wavelength = S(:,1);
    intensity = S(:,2);
    phase = S(:,3);
    if any(fftshift_indices == n) > 0
        intensities_without_mirrors(:,n) = fftshift(intensity);
        wavelengths_without_mirrors(:,n) = wavelength;
        phases_without_mirrors(:,n) = unwrap(-1.*fftshift(phase));
    else
        intensities_without_mirrors(:,n) = intensity;
        wavelengths_without_mirrors(:,n) = wavelength;
        phases_without_mirrors(:,n) = phase;
    end
    [H1(n)] = plot(AX(1),wavelengths_without_mirrors(:,n),intensities_without_mirrors(:,n));
    [H2(n)] = plot(AX(2),wavelengths_without_mirrors(:,n),phases_without_mirrors(:,n));
    set(H1(n),'Color',colors{n},'LineWidth',2)
    set(H2(n),'Color',colors{n},'LineStyle','--')
end
xlim(AX(1),[1000 1060])
xlim(AX(2),[1000 1060])
l = legend(H1,leg,'Location','South');
set(l,'FontSize',F_size_Label);
t = get(l,'Title');
set(t,'String','Measurement')

%##################### subplot(2,2,3) ####################
figure()
N = 2;

phase_diff = zeros(256,N);
for n=1:N
    phase_diff(:,n) = phases_without_mirrors(:,n)-phases_with_mirrors(:,n);
end
leg = {};
% CAREFUL !!! This is possibly incorrect since the wavelength values
% slightly vary. Subtracting a shifted from an unshifted array leads to
% small errors.
plot(wavelengths_without_mirrors(:,1),phase_diff(:,1),'r-','LineWidth',2)
leg{1} = sprintf(file_counter,1);
hold on
for n = 2:N
    plot(wavelengths_without_mirrors(:,n),phase_diff(:,n),'b-','LineWidth',2)
    leg{n} = sprintf(file_counter,n);
end
hold off
title('Phase difference $\varphi_{no mirror} - \varphi_{with mirror}$','FontSize',14,'Interpreter','LaTex')
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
l = legend(leg,'Location','South');
set(l,'FontSize',F_size_Label);
t = get(l,'Title');
set(t,'String','Measurement')

%##################### MIRROR SIGN CHECK ####################
% Calculate the Fourier transforms of the pulse with measured mirror phase
% added and subtracted. Compare the results with the retrieved pulse learn
% if the produced mirrors implement a positive or negative chirp.
% We know that the Menlo laser is positivle chirped, i.e. the phase curve
% will have a positive second derivative. Checking with the plots shows us
% that we have the correct phase for the Menlo retrieval.

file_end2 = '.bin.Ek.dat';
file_name = sprintf([file_folder file_counter without_mirrors file_end2],2);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
% The pulse does not have the correct phase and intensity yet. This comes
% from an ambiguity in the retrieval.
% In order to fix this we flip the time and multiply the phase by -1
% Time reversal corresponds to a sign flip of the phase --> Trebino p.128
time_without_mirrors = S(:,1);
temporal_intensity_without_mirrors = S(:,2);
temporal_phase_without_mirrors = S(:,3);
figure()
% subplot(2,2,1)
% [AX, H1, H2] = plotyy(time_without_mirrors,temporal_intensity_without_mirrors,time_without_mirrors,temporal_phase_without_mirrors);
% title(sprintf('The time flipped conjugated retrieved temporal pulse of the Menlo laser unchirped\nA negative quadratic coefficient in the temporal phase corresponds to a positive chirp'),'FontSize',14)
% set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
% set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
% set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)


% That's how the pulse should look like
wavelength_without_mirrors = wavelengths_without_mirrors(:,1);
spectral_intensity_without_mirrors = intensities_without_mirrors(:,1);
spectral_phase_without_mirrors = -1.*phases_without_mirrors(:,1);
[Int0,t0,Ek0] = compensation_calcFourierlimit(spectral_intensity_without_mirrors,wavelength_without_mirrors,spectral_phase_without_mirrors);
plotyy(t0,abs(Ek0).^2,t0,unwrap(angle(Ek0)))
title(sprintf('Fourier transform of the spectral pulse retrieval of the unchirped Menlo pulse\nData has been transformed to reflect the known physical properties\nof the unchirped Menlo output'),'FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('intesity [a.u.]','FontSize',F_size_Label)


figure()
file_end2 = '.bin.Ek.dat';
file_name = sprintf([file_folder file_counter with_mirrors file_end2],2);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
time_with_mirrors = S(:,1);
temporal_intensity_with_mirrors = S(:,2);
temporal_phase_with_mirrors = S(:,3);
[AX, H1, H2] = plotyy(time_with_mirrors,temporal_intensity_with_mirrors,time_with_mirrors,temporal_phase_with_mirrors);
title(sprintf('The retrieved temporal pulse of the chirped Menlo laser, after two bounces off the mirror'),'FontSize',14)
set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)


figure()
% We don't know if we measured the real compressor result pulse or the time
% flipped version of it. So we end up with four combinations we have to
% check in order to find out what the phase of our chirped mirror is:
%   Phi_mirror = Phi_comp - Phi_orig
%   Phi_mirror = -Phi_comp - Phi_orig
%   Phi_mirror = fliplr(-Phi_comp) - Phi_orig 
%   Phi_mirror = -fliplr(-Phi_comp) - Phi_orig
%--------------------------------------------------------------------------
Phi_mirror = interp1(Trubetskov_wavelength,Trubetskov_phase,wavelength_without_mirrors);
Phi_mirror(isnan(Phi_mirror)) = 0;
Phi_mirror = Phi_mirror .* 2; % two bounces
Phi_plus = spectral_phase_without_mirrors + Phi_mirror;
Phi_minus = spectral_phase_without_mirrors - Phi_mirror;

[Int_plus,t_plus,Ek_plus] = compensation_calcFourierlimit(spectral_intensity_without_mirrors,wavelength_without_mirrors,Phi_plus);
[Int_minus,t_minus,Ek_minus] = compensation_calcFourierlimit(spectral_intensity_without_mirrors,wavelength_without_mirrors,Phi_minus);
[Int_meas_plus,t_meas_plus,Ek_meas_plus] = compensation_calcFourierlimit(spectral_intensity_without_mirrors,wavelength_without_mirrors,spectral_phase_without_mirrors + (phase_diff(:,n)+line));
[Int_meas_minus,t_meas_minus,Ek_meas_minus] = compensation_calcFourierlimit(spectral_intensity_without_mirrors,wavelength_without_mirrors,spectral_phase_without_mirrors - (phase_diff(:,n)+line));
Int_FROG = abs(trapz(time_with_mirrors,temporal_intensity_with_mirrors));

H1 = zeros(1,5); % Allocate the array for the handles
H2 = zeros(1,5);
[AX,H1(1),H2(1)] = plotyy(t_plus,abs(Ek_plus).^2./Int_plus,t_plus,unwrap(angle(Ek_plus)));
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
[H1(2)] = plot(AX(1),t_minus,abs(Ek_minus).^2./Int_minus);
[H2(2)] = plot(AX(2),t_minus,unwrap(angle(Ek_minus)));
set(H1(2),'Color',colors{2},'LineWidth',2)
set(H2(2),'Color',colors{2},'LineStyle','--')
[H1(3)] = plot(AX(1),t_meas_plus,abs(Ek_meas_plus).^2./Int_meas_plus);
[H2(3)] = plot(AX(2),t_meas_plus,unwrap(angle(Ek_meas_plus)));
set(H1(3),'Color',colors{3},'LineWidth',2)
set(H2(3),'Color',colors{3},'LineStyle','--')
[H1(4)] = plot(AX(1),t_meas_minus,abs(Ek_meas_minus).^2./Int_meas_minus);
[H2(4)] = plot(AX(2),t_meas_minus,unwrap(angle(Ek_meas_minus)));
set(H1(4),'Color',colors{4},'LineWidth',2)
set(H2(4),'Color',colors{4},'LineStyle','--')
[H1(5)] = plot(AX(1),fliplr(time_with_mirrors')',temporal_intensity_with_mirrors./Int_FROG);
[H2(5)] = plot(AX(2),fliplr(time_with_mirrors')',-1.*temporal_phase_with_mirrors);
set(H1(5),'Color',colors{5},'LineWidth',2)
set(H2(5),'Color',colors{5},'LineStyle','--')
ylim(AX(1),[0 0.0013])
l = legend(H1,{'Phi0 + PhiTrubetskov' , 'Phi0 - PhiTrubetskov' , 'Phi0 + deltaPhiMeasured', 'Phi0 - deltaPhiMeasured' , 'FROG measurement'},'Location','NorthWest');
set(l,'FontSize',F_size_Label);


figure()
Phi_result = -phases_with_mirrors(:,1) - spectral_phase_without_mirrors; % + sign because the phase has to be multiplied by -1
plot(wavelength_without_mirrors,Phi_result,'b-')
title('Resulting spectral phase of the chirped mirror','FontSize',14)
xlabel('time [fs]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)


%##################### GD/GDD ####################
% We have to divide our measurement result by 2 because we had two bounces
% off the mirror
Phi_result = Phi_result ./ 2;
% A linear sum term within the phase has no physical meaning so we use this
% fact to make the calculated phase more alike the measured one
idx1_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1010);
idx1_M = auxiliary_find_closest_idx(wavelength_without_mirrors,1010);
idx2_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1050);
idx2_M = auxiliary_find_closest_idx(wavelength_without_mirrors,1050);
d1 = Trubetskov_phase(idx1_T) - Phi_result(idx1_M,1);
d2 = Trubetskov_phase(idx2_T) - Phi_result(idx2_M,1);
m = (d2-d1)/(1050-1010); % the gradient of our used line
line = m .* (wavelength_without_mirrors-1010) + Trubetskov_phase(idx1_T) - 17;
% Polynomial fit
poly_wavelength = wavelength_without_mirrors(idx2_M:idx1_M);
poly_phase = Phi_result(idx2_M:idx1_M,1);
poly_line = interp1(wavelength_without_mirrors,line,poly_wavelength);
poly = polyfit(poly_wavelength,poly_phase,6);
poly_fit = polyval(poly,wavelength_without_mirrors) + line;
poly = polyfit(wavelength_without_mirrors,poly_fit,6);
poly_fit = polyval(poly,wavelength_without_mirrors);
% We know the phase in wavelength domain -> transform to frequency space
% and calculate the GD and GDD
frequency_without_mirrors = (2*pi*c) ./ wavelength_without_mirrors; 
poly = polyfit(frequency_without_mirrors,poly_fit,6);
poly_deriv = polyder(poly);
poly_GD = polyval(poly_deriv,frequency_without_mirrors) .* -1e6; % ./1e9.*1e15 = .*1e6
poly_deriv = polyder(poly_deriv);
poly_GDD = polyval(poly_deriv,frequency_without_mirrors) .* -1e12; % ./1e18.*1e30 = .*1e12

figure() % New figure to show the results from the GD, GDD calculations
leg = {};
% CAREFUL !!! This is possibly incorrect since the wavelength values
% slightly vary. Subtracting a shifted from an unshifted array leads to
% small errors.
plot(wavelength_without_mirrors,Phi_result + line,'r-','LineWidth',2)
leg{1} = sprintf('Measured by T. Pleyer');
hold on
plot(Trubetskov_wavelength,Trubetskov_phase,'g-','LineWidth',2)
leg{2} = sprintf('Provided by M. Trubetskov');
plot(wavelength_without_mirrors,poly_fit,'b-','LineWidth',2)
leg{3} = sprintf('Polynomial fit to measured data');
hold off
xlim([1000 1060])
title('Measured and calculated phase curve of the produced chirped mirrors','FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
l = legend(leg,'Location','NorthWest');
set(l,'FontSize',F_size_Label);
% t = get(l,'Title');
% set(t,'String','Measurement')
figure()
% Before we can plot the data we have to make sure that they have the same
% offset from the axis in order to be comparable. The value of the
% theoretical measurement is taken as zero point
idx_end_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1050);
idx_end_V = auxiliary_find_closest_idx(Vova_wavelength,1050);
idx_end_M = auxiliary_find_closest_idx(wavelength_without_mirrors,1050);
d_V = Trubetskov_GD(idx_end_T) - Vova_GD(idx_end_V);
d_M = Trubetskov_GD(idx_end_T) - poly_GD(idx_end_M);
plot(Trubetskov_wavelength,Trubetskov_GD,'r-')
hold on
plot(Vova_wavelength,Vova_GD+d_V,'b-')
plot(wavelength_without_mirrors,poly_GD+d_M,'g-')
hold off
xlim([1010 1050])
title('Comparison of theoretical and measured GD values for the chirped mirrors','FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GD [fs]','FontSize',F_size_Label)
l = legend('Calculated by M. Trubetskov','Measured by E. Fedulova','Measured by T. Pleyer','Location','SouthEast');
set(l,'FontSize',F_size_Label);
figure()
plot(Trubetskov_wavelength,Trubetskov_GDD,'r-')
hold on
plot(Vova_wavelength,Vova_GDD,'b-')
plot(wavelength_without_mirrors,poly_GDD,'g-')
hold off
xlim([1010 1050])
title('Comparison of theoretical and measured GDD values for the chirped mirrors','FontSize',22)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GDD [fs^2]','FontSize',F_size_Label)
l = legend('Calculated by M. Trubetskov','Measured by E. Fedulova','Measured by T. Pleyer','Location','NorthWest');
set(l,'FontSize',F_size_Label);


% subplot(2,2,4)
