clear

dir_main = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/';
filename_base = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
endings       = {'center-flattened', ...
                 'only-third-bumb', ...
                 'smooth', ...
                 'wings-flattened', ...
                 'original', ...
                 'catmull', ...
                 'smooth-with-randn-pm05', ...
                 'smooth-with-randn-pm005', ...
                 'original-with-randn-pm05', ...
                 'original-with-randn-pm005', ...
                 'adjust-height-10241-10248', ...
                 'adjust-height-10268-10275'};

for i=1:length(endings)
% Daten\26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED\26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck-smooth.dat
inputfile = sprintf('%s%s.Speck-%s.dat',dir_main,filename_base,endings{i});

Sk = dlmread(inputfile);

Sk_cplx = sqrt(Sk(:,2)) .* exp(1i*Sk(:,3));
%Sk_cplx = Sk(:,4) - 1i*Sk(:,5);

[t2,Ek2] = Speck_Fourier(Sk(:,1).*1e-9,Sk_cplx); % Pass the arguments as SI units!

out(:,1) = t2*1e15; % Give the time as is in the Trebino format, i.e. fs
out(:,2) = abs(Ek2).^2;
out(:,2) = out(:,2)./max(out(:,2));
out(:,3) = unwrap(angle(Ek2));

dir_main2 = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/';

% Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Ek-smooth.dat
output_file = sprintf('%s%s.Ek-%s.dat',dir_main2,filename_base,endings{i});

dlmwrite(output_file,out,'delimiter','\t','precision',6)
end