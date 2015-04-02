%##################### INIT #####################
clear

folder_base = 'T:\LEX_measurements\hybrid data\';
folder_path = '20150112\Tobi\';

cd([folder_base folder_path]) % change to the working directory
styles = {'r-', 'b-', 'g-', 'c-', 'm-', 'k-'};
colors = {'red', 'blue', 'green', 'cyan', 'magenta', 'black'};
F_size_Label = 14;
c = 299792458;



%##################### SPECTRUM #####################
%Load the spectrum data from file into variable
initial_spec = dlmread('Measured\Spectrum_measured_Menlo_no_mirrors_1207022U1.TXT',';',8,0);
initial_spec = initial_spec(:,1:2);
resultant_spec = dlmread('Measured\Spectrum_measured_two_bounces_off_chirped_mirror_1207022U1.TXT',';',8,0);
resultant_spec = resultant_spec(:,1:2);

% Reduce the amount of unnecessary data (zeros)
lower_cutoff = 500;
upper_cutoff = 530;
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
subplot(2,2,1)
fftshift_indices = [1]; % These indicate which plots need an fftshift applied
leg = {}; % the legend object to be displayed
file_name = sprintf([file_folder file_counter with_mirrors file_end],1);
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
t = get(l,'Title');
set(t,'String','Measurement')
%legend([H1 H2],['Test1';'Test2'],['Muh';'Mäh'])



%##################### subplot(2,2,2) ####################
N = 2; % Number of measurements made
subplot(2,2,2)
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
t = get(l,'Title');
set(t,'String','Measurement')

%##################### subplot(2,2,3) ####################
subplot(2,2,3)
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
t = get(l,'Title');
set(t,'String','Measurement')



%##################### GD/GDD ####################
% We have to divide our measurement result by 2 because we had two bounces
% off the mirror
phase_diff = phase_diff ./ 2;
% Get the theoretical and measured values for the GD and GDD
Trubetskov = dlmread('Trubetskov_Mirror_Specification.txt','\t',1,0);
Trubetskov_wavelength = Trubetskov(:,1);
Trubetskov_phase = Trubetskov(:,2);
% A linear sum term within the phase has no physical meaning so we use this
% fact to make the calculated phase more alike the measured one
idx1_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1010);
idx1_M = auxiliary_find_closest_idx(wavelengths_with_mirrors(:,1),1010);
idx2_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1050);
idx2_M = auxiliary_find_closest_idx(wavelengths_with_mirrors(:,1),1050);
d1 = Trubetskov_phase(idx1_T) - phase_diff(idx1_M,1);
d2 = Trubetskov_phase(idx2_T) - phase_diff(idx2_M,1);
m = (d2-d1)/(1050-1010); % the gradient of our used line
line = m .* (wavelengths_with_mirrors(:,1)-1010) + Trubetskov_phase(idx1_T) - 17;
Trubetskov_GD = Trubetskov(:,3);
Trubetskov_GDD = Trubetskov(:,4);
Vova = dlmread('Vova_Mirror_Data_Narrow_Band.dat','\t',2,0);
Vova_wavelength = Vova(:,1);
Vova_GD = Vova(:,2);
Vova_GDD = Vova(:,3);

% Calculate the values for GD and GDD as I measured them
GD = zeros(255,N);
GDD = zeros(254,N);
for n=1:N
    w = wavelengths_without_mirrors(:,n) .* 1e-9;
    p = phase_diff(:,n);
    GD(:,n) = diff(p)./(-diff(w));
    GDD(:,n) = diff(GD(:,n))./(-diff(w(1:255)));
    GD(:,n) = (-1) .* (w(1:255)).^2 ./ (2*pi*c) .* GD(:,n) .* 1e15;
    GD(:,n) = smooth(GD(:,n));
    GDD(:,n) = (w(1:254)).^3 ./ (2*pi*c.^2) .* GDD(:,n) .* 1e30;
    GDD(:,n) = smooth(GDD(:,n));
end

figure() % New figure to show the results from the GD, GDD calculations
subplot(2,2,1)
leg = {};
% CAREFUL !!! This is possibly incorrect since the wavelength values
% slightly vary. Subtracting a shifted from an unshifted array leads to
% small errors.
plot(wavelengths_with_mirrors(:,1),phase_diff(:,1)+line,'r-','LineWidth',2)
legend_string = ['Measurement ' file_counter ' with an additional\nbias of the form a+b*lambda'];
leg{1} = sprintf(legend_string,1);
hold on
for n = 2:N
    plot(wavelengths_with_mirrors(:,n),phase_diff(:,n)+line,'b-','LineWidth',2)
    leg{n} = sprintf(legend_string,n);
end
plot(Trubetskov_wavelength,Trubetskov_phase,'g-','LineWidth',2)
leg{N+1} = sprintf('Provided by M. Trubetskov');
hold off
xlim([1000 1060])
title('Phase difference between direct Menlo Output and Chirped Mirror compression','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('phase [rad]','FontSize',F_size_Label)
l = legend(leg,'Location','NorthWest');
% t = get(l,'Title');
% set(t,'String','Measurement')
subplot(2,2,2)
% Before we can plot the data we have to make sure that they have the same
% offset from the axis in order to be comparable. The value of the
% theoretical measurement is taken as zero point
idx_end_T = auxiliary_find_closest_idx(Trubetskov_wavelength,1050);
idx_end_V = auxiliary_find_closest_idx(Vova_wavelength,1050);
idx_end_M = auxiliary_find_closest_idx(wavelengths_with_mirrors(:,1),1050);
d_V = Trubetskov_GD(idx_end_T) - Vova_GD(idx_end_V);
d_M = Trubetskov_GD(idx_end_T) - GD(idx_end_M,1);
plot(Trubetskov_wavelength,Trubetskov_GD,'r-')
hold on
plot(Vova_wavelength,Vova_GD+d_V,'b-')
plot(wavelengths_with_mirrors(1:255,1),GD(:,1)+d_M,'g-')
hold off
xlim([1010 1050])
title('Comparison of theoretical and measured GD values for the chirped mirrors','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GD [fs]','FontSize',F_size_Label)
legend('Calculated by M. Trubetskov','Measured by E. Fedulova','Measured by T. Pleyer','Location','SouthEast')
subplot(2,2,3)
plot(Trubetskov_wavelength,Trubetskov_GDD,'r-')
hold on
plot(Vova_wavelength,Vova_GDD,'b-')
hold off
xlim([1010 1050])
title('Comparison of theoretical and measured GDD values for the chirped mirrors','FontSize',14)
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('GDD [fs^2]','FontSize',F_size_Label)
legend('Calculated by M. Trubetskov','Measured by E. Fedulova','Measured by T. Pleyer','Location','NorthWest')
subplot(2,2,4)
% Find the GDD via a polynomial fit and analytical derivative
poly_wavelength = wavelengths_with_mirrors(idx2_M:idx1_M,1);
poly_phase = phase_diff(idx2_M:idx1_M,1);
poly_line = m .* (poly_wavelength-1010) + Trubetskov_phase(idx1_T) - 17;
poly = polyfit(poly_wavelength,poly_phase,6);
poly_deriv = polyder(polyder(poly));
poly_GDD = polyval(poly_deriv,poly_wavelength) .* (poly_wavelength .* 1e-9).^3 ./ (2*pi*c.^2) .* 1e30;