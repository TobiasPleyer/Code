%%Thorlabs grating high efficiency
%%GR50-0310 300 Groves/mm angle 8°36'
clear

%% Init
GratingName = 'test';%'Rchdsn_53 590R';
N = 300;%grove density #/mm
alpha = 15;%8.6;
alpha = 2*pi/360*(alpha); %angle of incidence

filename_base = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
dir_main = '\\gar-sv-home01\Tobias.Pleyer\Desktop\Matlab Moritz\FROG\';

range_time = '[-5000:5000]';
range_spec = '[1010:1050]';

N_pictures = 1;
D_min = -20e-2;                                                            % The minimum distance between the two gratings. A negative value corresponds to a virtual negative distance as it could be introduced by a spectrometer
D_max = 11e-2;                                                             % The minimum distance between the two gratings.
phase_factor = -1;
NPhaseFit = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_spec = sprintf('%sSpectrum\',dir_main);
dir_time = sprintf('%sTime\',dir_main);
dir_data = sprintf('%sData\',dir_main);
Npoints = 100;
%%
    c = 299792458; %m/s

Dist = linspace(D_min,D_max, N_pictures);                                  % This creates the range of grating distances that are run through

for count_main = 1:N_pictures
    
    filenameextra = sprintf('%03d_',count_main);
    
    D0 = Dist(count_main); 

    D = linspace(0,5.*D0,Npoints+1);                                       % distance of gratings
    ND0 = Npoints/5+1;

    lambda = linspace(1020e-9,1040e-9,400);
    lambda_0 = 1030e-9;                                                    % central wavelength in m
    omega = 2.*pi.*c./lambda;
    omega_0 = 2.*pi.*c/lambda_0;
    beta_0 = -asin(2.*pi.*c.*N.*1e3./omega_0 - sin(alpha));

    D1 = 2.*D/c;
    D2 = -2.*D./c.*(lambda_0./(2.*pi.*c)).* ...
         (lambda_0.*N.*1e3./(cos(beta_0))).^2;
    D3 = -3.*D2.*lambda_0./(2.*pi.*c).* ...
         (1+lambda_0.*N.*1e3.*sin(beta_0)./(cos(beta_0)).^2);
    D4 = 3.*D2.*(lambda_0./(2.*pi.*c)).^2.* ...
         (4.*(1+lambda_0.*N.*1e3.*sin(beta_0)./(cos(beta_0)).^2).^ ...
         2+(lambda_0.*N.*1e3/(cos(beta_0)).^2).^2);

    %plot(D*1e2,D2.*1e30)

    beta_omega = -asin(2.*pi.*c.*N.*1e3./omega - sin(alpha));

    phi = 2.*omega./c.*D0.*cos(beta_0-beta_omega)-2.*omega_0./c.* ...
          D(ND0) - D1(ND0).*(omega-omega_0);
    % The Taylor series of the phase up to a specific order: 
    % Phi = Sum 1/k! d^kPhi/domega^k * (omega-omega0)^k
    phi_tay = 0.*2.*omega_0./c.*D(ND0)...                                  % constant frequencyl. Not of interest for us
             + 0.*D1(ND0).*(omega-omega_0)...                              % Group Delay. Not of interest for us. Just causes a temporal shift 
             + 1/2.*D2(ND0).*(omega-omega_0).^2 ...                        % GDD
             + 1/6.*D3(ND0).*(omega-omega_0).^3 ...                        % TOD
             + 1/24 .*D4(ND0).*(omega-omega_0).^4;                         % FOD
    plot(omega, phi, omega, phi_tay)


    %%Apply to data
    
    filename_Speck = sprintf('%s%s.Speck.dat',dir_main, filename_base);
    Speck = dlmread(filename_Speck);

    filename_Ek = sprintf('%s%s.Ek.dat',dir_main, filename_base);
    Ek_orig = dlmread(filename_Ek);
    Ek_orig = Ek_orig(:,1:3);

    f = c./(Speck(:,1).*1e-9);
    omega = 2.*pi.*f;
    T = 1/(f(2)-f(1));
    t = -T/2.*linspace(-1,1,size(f,1));
    t=t';

    beta_omega = -asin(2.*pi.*c.*N.*1e3./omega - sin(alpha));
    phi_grating = 2.*omega./c.*D0.*cos(beta_0-beta_omega)-2.*omega_0./ ...
                  c.*D(ND0)- D1(ND0).*(omega-omega_0);

    p = polyfit((omega-omega_0).*1e-15, phi_grating, NPhaseFit);
    phifit = polyval(p, (omega-omega_0).*1e-15);
    
    Sw = (Speck(:,1).*1e-9).^2./(2*pi*c).*Speck(:,2);
    Sw = Sw./max(Sw);

    Phiw = Speck(:,3).*phase_factor;

    A_1f = sqrt(Sw).*exp(-1i.*Phiw);

    H_ef = exp(-1i.*phi_grating);
    A_2f = A_1f.*H_ef;

    Ek=fftshift(fft(ifftshift(A_2f)));
    I = abs(Ek/max(abs(Ek))).^2;

    Ek_compr = Ek;
    I_compr = I;
    Speck_compr = A_2f;
    FWHM_pulse = abs(fwhm_moritz(t.*1e15,I_compr)*1e-3);


    %%
    %Write data to file and plot
    filen_EkGrating = sprintf('%s%s.Ek_%s_D=%05.1fcm.dat',dir_data,filename_base,GratingName,D0.*1e2);
    dlmwrite(filen_EkGrating,[t*1e15, I_compr, -unwrap(angle(Ek_compr)), ...
              real(Ek_compr), imag(Ek_compr)],'\t');

    filen = sprintf('%s.Ek_%s_D=%05.1fcm',filename_base,GratingName,D0.*1e2);
    % %write gnuplot file (reconstructed)
    gnu_file_plt=sprintf('%s%s.plt',dir_data,filen);
    fid_gnu = fopen(gnu_file_plt,'wt');
    fprintf(fid_gnu,['set term jpeg\n'...
                    'set output ''%s%s%s.jpeg''\n'...
                    'set title ''Time Domain FWHM=%.3fps N=%03.u#/mm Alpha_in=%.1f°''\n'...
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
                    dir_time,filenameextra,filen,FWHM_pulse,N, 360/2/pi*alpha,range_time ,filen_EkGrating,filen_EkGrating);
    fclose(fid_gnu);

    %execute gnuplot to create jpeg
    gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
    eval(['!gnuplot.exe', gnu_file_plt_eval]);

    %%
    %Write data to file and plot
    filen_SpeckGrating = sprintf('%s%s.Speck%s_D=%05.1fcm.dat',dir_data,filename_base,GratingName,D0.*1e2);
    dlmwrite(filen_SpeckGrating,[c./f.*1e9, ...
        abs(Speck_compr).^2*2*pi.*f.^2./c/max(abs(Speck_compr).^2*2*pi.*f.^2./c),...
        -unwrap(angle(Speck_compr)), real(Speck_compr), imag(Speck_compr)],'\t');

    filen = sprintf('%s.Speck%s_D=%05.1fcm',filename_base,GratingName,D0.*1e2);
    % %write gnuplot file (reconstructed)
    gnu_file_plt=sprintf('%s%s.plt',dir_data,filen);
    fid_gnu = fopen(gnu_file_plt,'wt');
    fprintf(fid_gnu,['set term jpeg\n'...
                    'set output ''%s%s%s.jpeg''\n'...
                    'set title ''Spectral Domain FWHM=%.3fps N=%03.u#/mm Alpha_in=%.1f°''\n'...
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
                    dir_spec,filenameextra, filen,FWHM_pulse, N, 360/2/pi*alpha,range_spec, filen_SpeckGrating,filen_SpeckGrating);
    fclose(fid_gnu);

    %execute gnuplot to create jpeg
    gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
    eval(['!"gnuplot.exe"', gnu_file_plt_eval]);
    
    [ax h1 h2] = plotyy(omega, Sw, omega, Phiw);
    hold(ax(2), 'on');
    plot(ax(2), omega, -phifit);
    plot(ax(2), omega, Phiw+phifit,'r');
    hold(ax(2), 'off')

end

sprintf('Dispersion of grating D=%.1fcm\n',D0*100)
sprintf('%02d: %.1e\n',[NPhaseFit:-1:0; p])
