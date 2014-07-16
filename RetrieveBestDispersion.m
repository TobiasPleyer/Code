%%Script reads a FROG file and tries to find GDD, TOD and FOD needed
% for compression. The respective GDD, TOD and FOD are obtained as the
% coefficients of a simple polynomial fit to the measured phase of the
% pulse.


clear

%% Inititialize
debug=false;
dir_main = '\\gar-sv-home01\Tobias.Pleyer\Desktop\Matlab Moritz\FROG';
dir_spec = 'FROG';
dir_time = 'Time/';
dir_data = 'Data/';
filename_base = '04_AIR_FROG_30.4W_RTT=7.415us_Ip=12.1A.bin';
%filename_base = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin_D=025.0cm';
%filename_base = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';

ratioSet = 5;

manual = false;
pman = [-1e7 1.3e6 -5e4 1e2 -5.6];
NPhaseFit = 5;
FitUseOnly = NPhaseFit;

Threshold = 0.005;

correct_phase_unwrap = [137, 0];

lambda0 = 1030e-9;                                                         % central wavelength in m
lambdaMin = 1020e-9;
lambdaMax = 1040e-9;
Nlambda = 400;

phase_factor = 1;
range_time = '[-5000:5000]';
range_spec = '[1010:1050]';
filename_extra = sprintf('Order=%01d_Man=%d_',NPhaseFit,manual);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458; %m/s

dir_spec = sprintf('%s%s',dir_main, dir_spec);
dir_time = sprintf('%s%s',dir_main, dir_time);
dir_data = sprintf('%s%s',dir_main, dir_data);


%lambda = linspace(lambdaMin,lambdaMax,Nlambda);
omega0 = 2.*pi.*c/lambda0;

%% read data

filename_Speck = sprintf('%s%s.Speck.dat',dir_main, filename_base);
Speck = dlmread(filename_Speck);                                           % read the csv file created by the FROG algorithm
f = c./(Speck(:,1).*1e-9);                                                 % the frequencies corresponding to the given wavelengths
omega = 2.*pi.*f;                                                          % the angular frequencies corresponding to the respective frequencies
T = 1/(f(2)-f(1));                                                         % find the scale in the time domain that spans the whole pulse (borders)
t = T/2.*linspace(-1,1,size(f,1));                                         % create the time domain range
t=t';

Sw = (Speck(:,1).*1e-9).^2./(2*pi*c).*Speck(:,2);                          % Find the spectral intensity distribution and adds a frequency dependent factor in order to maintain equal integrals
Sw = Sw./max(Sw);                                                          % Normalize the spectrum
Phiw = Speck(:,3);                                                         % The phase factor of the respective frequency
Phiw = Phiw.*phase_factor;
for i = 1:size(correct_phase_unwrap,1)                                     % If a look at the graph reveals some sort phase offset from a certain
    ind = correct_phase_unwrap(i,1);                                       % frequency component on, it is possible to correct for that by 
    Phiw =  [Phiw(1:ind); Phiw(ind+1:end)+correct_phase_unwrap(i,2)*2*pi]; % changing the correct_phase_unwrap[freq,multiple of 2pi] variable.
    if (debug)                                                             % The parameter i defaults to 1, but can be >1 if correct_phase_unwrap
        plot(omega,Phiw)                                                   % is extended by more frequency adjustments, e.g.
        hold on                                                            % correct_phase_unwrap = [137, 0; 150,1];
        plot(omega,Speck(:,3));                                            % The second parameter always constitutes the phase change to be applied
        hold off                                                           % which will be applied to all frequency components following the one
        k = waitforbuttonpress;                                            % provided as the first argument.
    end
end

A_1f = sqrt(Sw).*exp(-1i.*Phiw);                                           % This constitutes the complex electric field in the spectral domain
%Get Fourier limit. The Fourier limit of a pulse equals the pulses
%intensity with a flat, i.e. constant, phase. In our case phi=0.
EkFourier=fftshift(fft(ifftshift(sqrt(Sw))));                              % Find the temporal electric field via fourier transform. The fftshifts are applied in order to keep the data symmetrically centered
IFourier = abs(EkFourier/max(abs(EkFourier))).^2;                          % Find the intensity of the electric field
FWHMFourier = abs(fwhm_moritz(t.*1e15,IFourier)*1e-3);                     % Find the full width half maximum duration of the obtained temporal pulse (may be an inappropriate measure)
RMSFourier = RMS_moritz(t.*1e15,IFourier)*1e-3;                            % Find the root mean square duration of the obtained temporal pulse (may be an inappropriate measure)
EquivFourier = Equiv_moritz(t.*1e15,IFourier)*1e-3;                        % Find equivalent fourier limit duration of the obtained temporal pulse (may be an inappropriate measure)


filename_Ek = sprintf('%s%s.Ek.dat',dir_main, filename_base);
Ek_orig = dlmread(filename_Ek);
Ek_orig = Ek_orig(:,1:3);




%Ek=fftshift(fft(ifftshift(A_2f)));
%I = abs(Ek/max(abs(Ek))).^2;

IndRelev(1) = find(Sw > Threshold, 1, 'first');                            % This defines the lower limit of the range of interest for us, where the intensity is !=0
IndRelev(2) = find(Sw > Threshold, 1, 'last');                             % This defines the upper limit of the range of interest for us, where the intensity is !=0

PhiFit = Phiw(IndRelev(1):IndRelev(2));                                    % Cut out the phase range we want to consider from the original phase according to our above limits

omega_fit = (omega(IndRelev(1):IndRelev(2))-omega0).*1e-15;                % Do the same for the corresponding angular frequency range and give it relativ to the central angular frequency

if (manual)
    p = pman;
else
    p = polyfit(omega_fit, PhiFit, NPhaseFit);                             % Find a polynomial fit of the desired degree NPhaseFit
end


phi_fit = 0;
phi_add = 0;
for i = NPhaseFit-FitUseOnly+1:NPhaseFit+1                                 % Create the fitted polynomial for the present data
    phi_fit = phi_fit+p(i).*omega_fit.^(NPhaseFit-(i-1));                  % This is the phase we will be plotting for visual comparison
    phi_add = phi_add+p(i).*((omega-omega0).*1e-15).^(NPhaseFit-(i-1));    % This is the phase term that will be added to every phase to counteract the present phase in an attempt to flatten the latter
end


H_ef = exp(1i.*phi_add);                                                   % This is the term that will be multiplied to the existing electric field to transform its phase
A_2f = A_1f.*H_ef;                                                         % A_1f.*H_ef = sqrt(Sw).*exp(-1i.*Phiw).*exp(1i.*phi_add) = sqrt(Sw).*exp(-1i.*(Phiw1i-phi_add)) ~ sqrt(Sw).*exp(constant)

Ek=fftshift(fft(ifftshift(A_2f)));                                         % Create the Fourier transform of the field with the compensated phase
I = abs(Ek/max(abs(Ek))).^2;                                               % Find the intensity of the compensated field

FWHMPulse = abs(fwhm_moritz(t.*1e15,I)*1e-3);                              % Find the full width half maximum duration of the obtained temporal pulse (may be an inappropriate measure)

%% Write data to file and plot
filen_Ek = sprintf('%s%s%s.Ek.dat',...
    dir_data,filename_extra,filename_base);
dlmwrite(filen_Ek,[t*1e15, I, -unwrap(angle(Ek)), ...
                real(Ek), imag(Ek)],'\t');

filen = sprintf('%s%s.Ek',filename_extra,filename_base);
% %write gnuplot file (reconstructed)
gnu_file_plt=sprintf('%s%s.plt',dir_data,filen);
fid_gnu = fopen(gnu_file_plt,'wt');
fprintf(fid_gnu,['set term jpeg\n'...
    'set output ''%s%s.jpeg''\n'...
    'set title ''Time Domain FWHM=%.3fps''\n'...
    'set style line 1 lt 1 lw 1\n'... 
    'set style line 2 lt 2 lw 1\n'...
    'set y2label ''Phase in [rad]''\n'...
    'set ytics nomirror\n'...
    'set y2tics auto\n'...
    'set xtics auto\n'...
    'set yrange [-0.05:1.05]\n'...
    'set y2range [:]\n'...
    'set ytics 0,1.0/4\n'...
    'set xlabel ''Time [fs]''\n'...
    'set ylabel ''Intensity (a.u)''\n'...
    'set xrange %s\n'...
    'plot ''%s'' using 1:2 axis x1y1 title ''Envelope'' with lines ls 1,\\\n'...
    '''%s'' using 1:3 axis x1y2 title ''Phase'' with lines ls 2\n'...
    'set output '],...
    dir_time,filen,FWHMPulse,range_time ,filen_Ek,filen_Ek);
fclose(fid_gnu);

%execute gnuplot to create jpeg
gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
eval(['!gnuplot.exe', gnu_file_plt_eval]);

%%
%Write data to file and plot
filen_Speck = sprintf('%s%s%s.Speck.dat',...
    dir_data,filename_extra,filename_base);
 dlmwrite(filen_Speck,[c./f.*1e9, ...
            abs(A_2f).^2*2*pi.*f.^2./c/max(abs(A_2f).^2*2*pi.*f.^2./c),...
            -unwrap(angle(A_2f)), real(A_2f), imag(A_2f)],'\t');

filen = sprintf('%s%s.Speck',filename_extra,filename_base);
% %write gnuplot file (reconstructed)
gnu_file_plt=sprintf('%s%s.plt',dir_data,filen);
fid_gnu = fopen(gnu_file_plt,'wt');
fprintf(fid_gnu,['set term jpeg\n'...
    'set output ''%s%s.jpeg''\n'...
    'set title ''Spectral Domain FWHM=%.3fps''\n'...
    'set style line 1 lt 1 lw 1\n'... 
    'set style line 2 lt 2 lw 1\n'...
    'set y2label ''Phase in [rad]''\n'...
    'set ytics nomirror\n'...
    'set y2tics auto\n'...
    'set xtics auto\n'...
    'set yrange [-0.05:1.05]\n'...
    'set y2range [:]\n'...
    'set ytics 0,1.0/4\n'...
    'set xlabel ''Wavelength [nm]''\n'...
    'set ylabel ''Intensity (a.u)''\n'...
    'set xrange %s\n'...
    'plot ''%s'' using 1:2 axis x1y1 title ''Envelope'' with lines ls 1,\\\n'...
    '''%s'' using 1:3 axis x1y2 title ''Phase'' with lines ls 2\n'...
    'set output '],...
    dir_spec,filen,FWHMPulse,range_spec, filen_Speck,filen_Speck);
fclose(fid_gnu);

%execute gnuplot to create jpeg
gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
eval(['!"gnuplot.exe"', gnu_file_plt_eval]);

figure(1)
plot(t*1e15, I, 'b','LineWidth',2)
hold on
plot(t*1e15, IFourier, 'r')
plot(t*1e15, Ek_orig(:,2), 'g')
hold off
xlabel('Time [fs]');
ylabel('Intensity (a.u.)');

figure(2)
[ax, h1, h2] = plotyy(omega, Sw, omega, Phiw);
hold(ax(2),'on')
plot(ax(2), omega_fit.*1e15+omega0, phi_fit, 'r');
plot(ax(2), omega, Phiw-phi_add, 'g');
hold(ax(2),'off')
xlabel('Frequency [Hz]');
ylabel(ax(1),'Intensity (a.u.)');
ylabel(ax(2),'Phase (rad)');

readme_file = sprintf('%s%s.txt',dir_data,filen);
fid_readme = fopen(readme_file,'wt');
fprintf(fid_readme, 'Ordes of fit fs^n\n');

N=length(p);
for i=0:N-1
    fitorder(N-i) = p(N-i)*factorial(i);
end
fprintf(fid_readme,'%02d: %.1e\n', [NPhaseFit:-1:0; fitorder]);
fclose(fid_readme);