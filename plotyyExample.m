%##################### INIT #####################
folder_base = 'T:\LEX_measurements\hybrid data\';
folder_path = '20150112\Tobi\';

cd([folder_base folder_path]) % change to the working directory
styles = {'r-', 'b-', 'g-', 'c-', 'm-', 'k-'};
colors = {'red', 'blue', 'green', 'cyan', 'magenta', 'black'};
F_size_Label = 14;


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
subplot(2,2,1)
plot(initial_spec(:,1),initial_spec(:,2)./max(initial_spec(:,2)),'b-')
hold on
plot(resultant_spec(:,1),resultant_spec(:,2)./max(resultant_spec(:,2)),'r-')
plot([507.5 507.5 527.5 527.5],[0 1 1 0],'g--')
ylim([0,1.1])
xlabel('wavelength [nm]','FontSize',F_size_Label)
ylabel('Intensity [$\frac{N_{samples}}{max(N_{samples})}$]','Interpreter','Latex','FontSize',F_size_Label+6)
title(sprintf('Comparison between the initial spectrum of the Menlo laser and\nthe spectrum after two bounces off the chirped mirror (after SHG crystal)'),'FontSize',14)
legend('initial spectrum','after two bounces off mirror','design limits of mirror','Location','NorthWest')
hold off


%##################### FROG #####################
file_base1 = '_SHG-FROG_Compressor_angle_check_Reference_M1';
file_base2 = '_SHG-FROG_Menlo_no_mirrors';
file_counter = '%03d';
file_end = '.bin.Speck.dat';
file_folder = 'FROG/';

%##################### subplot(2,2,2) #####################
N = 2; % Number of measurements made
subplot(2,2,2)
fftshift_indices = [1]; % These indicate which plots need an fftshift applied
leg = {}; % the legend object to be displayed
file_name = sprintf([file_folder file_counter file_base1 file_end],1);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
leg{1} = sprintf(file_counter,1);
wavelength = S(:,1)/2;
intensity = S(:,2);
phase = S(:,3);
H1 = zeros(1,N); % Allocate the array for the handles
H2 = zeros(1,N);
if any(fftshift_indices == 1) > 0
    [AX,H1(1),H2(1)] = plotyy(wavelength,fftshift(intensity),wavelength,unwrap(-1.*fftshift(phase)));
else
    [AX,H1(1),H2(1)] = plotyy(wavelength,intensity,wavelength,phase);
end
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
    file_name = sprintf([file_folder file_counter file_base1 file_end],n);
    disp(file_name)
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    leg{n} = sprintf(file_counter,n);
    wavelength = S(:,1)/2;
    intensity = S(:,2);
    phase = S(:,3);
    if any(fftshift_indices == n) > 0
        [H1(n)] = plot(AX(1),wavelength,fftshift(intensity));
        [H2(n)] = plot(AX(2),wavelength,unwrap(-1.*fftshift(phase)));
    else
        [H1(n)] = plot(AX(1),wavelength,intensity);
        [H2(n)] = plot(AX(2),wavelength,phase);
    end
    set(H1(n),'Color',colors{n},'LineWidth',2)
    set(H2(n),'Color',colors{n},'LineStyle','--')
end
xlim(AX(1),[500 530])
xlim(AX(2),[500 530])
l = legend(H1,leg,'Location','South');
t = get(l,'Title');
set(t,'String','Measurement')
%legend([H1 H2],['Test1';'Test2'],['Muh';'Mäh'])

%##################### subplot(2,2,3) ####################
N = 2; % Number of measurements made
subplot(2,2,3)
fftshift_indices = [1 2]; % These indicate which plots need an fftshift applied
%leg = {}; % the legend object to be displayed
file_name = sprintf([file_folder file_counter file_base2 file_end],1);
disp(file_name)
try
    S = dlmread(file_name);
catch err
    fprintf('Error: File not found\n')
end
leg{1} = sprintf(file_counter,1);
wavelength = S(:,1)/2;
intensity = S(:,2);
phase = S(:,3);
H1 = zeros(1,N); % Allocate the array for the handles
H2 = zeros(1,N);
if any(fftshift_indices == 1) > 0
    [AX,H1(1),H2(1)] = plotyy(wavelength,fftshift(intensity),wavelength,unwrap(-1.*fftshift(phase)));
else
    [AX,H1(1),H2(1)] = plotyy(wavelength,intensity,wavelength,phase);
end
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
    file_name = sprintf([file_folder file_counter file_base2 file_end],n);
    disp(file_name)
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    leg{n} = sprintf(file_counter,n);
    wavelength = S(:,1)/2;
    intensity = S(:,2);
    phase = S(:,3);
    if any(fftshift_indices == n) > 0
        [H1(n)] = plot(AX(1),wavelength,fftshift(intensity));
        [H2(n)] = plot(AX(2),wavelength,unwrap(-1.*fftshift(phase)));
    else
        [H1(n)] = plot(AX(1),wavelength,intensity);
        [H2(n)] = plot(AX(2),wavelength,phase);
    end
    set(H1(n),'Color',colors{n},'LineWidth',2)
    set(H2(n),'Color',colors{n},'LineStyle','--')
end
xlim(AX(1),[500 530])
xlim(AX(2),[500 530])
l = legend(H1,leg,'Location','South');
t = get(l,'Title');
set(t,'String','Measurement')

%The y-scale is calculated from all availbale x-data, not just the one
%used by xlim -> find the right scaling
% idx1 = auxiliary_find_closest_idx(t2,-1500);
% idx2 = auxiliary_find_closest_idx(t2,1500);
% p =unwrap(-angle(Ek2));
% m = floor(min(p(min(idx1,idx2):max(idx1,idx2))));
% M = ceil(max(p(min(idx1,idx2):max(idx1,idx2))));
% ylim(AX(2),[m M])
% set(AX(2),'YTick',m:5:M)
% set(AX(1),'Box','off')