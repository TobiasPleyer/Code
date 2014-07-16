clear;

dir_main      = 'Daten/';
data_base     = '26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/';
%picture_base  = 'Pictures\26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_PICTURES\';
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
                 'adjust-height-10268-10275', ...
                 'trebino'};

for i=1:length(endings)
filename_Ek = sprintf('%s%s%s.Ek-%s.dat',dir_main,data_base,filename_base,endings{i});
Ek_orig = dlmread(filename_Ek);
Ek = Ek_orig(:,1:3);

lambda0 = 1030;
l_min = 990;
l_max = 1060;
gridSize = 1000;
tau_range = 7000;
tbuf = 000;

%[t, f, spec, type] = SHG_FROG(Ek, tau_range, tbuf, lambda0, l_min, l_max, gridSize);
[t, f, spec, type] = PG_FROG(Ek, tau_range, tbuf, lambda0, l_min, l_max, gridSize);

[TAU, F] = meshgrid(t, f);


% figure(1)
% surf(TAU, F, spec,'EdgeColor','none');
% view(0,90)
% title(type)
% xlabel('Time [fs]')
% ylabel('Frequency [1/fs]')

figure_name = sprintf('%s-%s','26-AIR-FROG-31W-RTT=7us-Ip=12A-Ek',endings{i});
figure_file = sprintf('%s-%s','Daten/Pictures/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_PICTURES/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek',endings{i});
figure(2)

[ax,h1,h2] = plotyy(Ek(:,1),Ek(:,2),Ek(:,1),Ek(:,3));
legend([h1,h2],'intensity','phase')
title(figure_name)
saveas(gcf,figure_file,'jpg')

gnu_file_dat_1 = sprintf('%s-%s.A',filename_base,endings{i});
fullfilename = sprintf('%sDat/%s_%s_gnuplot.dat',dir_main,type,gnu_file_dat_1);
savegpbin_moritz(t, f , spec', fullfilename);

%write gnuplot file (measured)
gnu_file_plt=sprintf('%sFrog/%s_%s.plt',dir_main,type,gnu_file_dat_1);
fid_gnu_1 = fopen(gnu_file_plt,'wt');
fprintf(fid_gnu_1,['set term jpeg\n'...
        'set output ''%sPictures/%s/%s.jpeg''\n'],dir_main,type,gnu_file_dat_1);
str = sprintf(['set pm3d map\n'...
    'set palette defined (0 "black", 0.01 "red",0.16 "yellow",0.46 "green", 0.81 "blue", 1 "white")\n'...
    'set title ''%s trace for %s\n'...
    'set ytics auto\n'...
    'set cblabel ''Intensity (a.u)''\n'...
    'set xtics auto\n'...
    'set yrange [:]\n'...
    'set xlabel ''Delay [fs]''\n'...
    'set ylabel ''Wavelength [nm]''\n'...
    'set xrange [:]\n'...
    'splot ''%s'' binary matrix with pm3d notitle\n'],...
    type,figure_name, fullfilename);
fprintf(fid_gnu_1,'%s',str);
fprintf(fid_gnu_1,'set output ');

fclose(fid_gnu_1);

%execute gnuplot to create jpeg
gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
eval(['!"C:/Program Files (x86)/gnuplot/bin/gnuplot.exe"', gnu_file_plt_eval]);
end


% filename_recon = sprintf('%s%s%s.Ek-%s.dat',dir_main,data_base,filename_base,endings{6});
% filename_treb  = sprintf('%s%s%s.Ek-%s.dat',dir_main,data_base,filename_base,endings{7});
% Ek_recon = dlmread(filename_recon);
% Ek_treb  = dlmread(filename_treb);
% 
% wavelength1 = Ek_recon(:,1);
% wavelength2 = Ek_treb(:,1); %reshape(fliplr(reshape(Ek_treb(:,1),1,256)),256,1);
% scaling_factor = 1; %0.809226/0.93;
% I_recon = Ek_recon(:,2)*scaling_factor;
% P_recon = Ek_recon(:,3);
% I_treb = Ek_treb(:,2);
% P_treb = Ek_treb(:,3);
% 
% figure(2)
% [ax,h1,h2] = plotyy([wavelength1,wavelength2],[I_recon I_treb],[wavelength1,wavelength2],[P_recon P_treb]);
% text(2000,0.85,sprintf('recon scaled by %f',scaling_factor))
% legend([h1;h2],'reconstructed intensity','Trebino intensity','reconstructed phase','Trebino phase','Location','NorthWest')
% title('Comparison between the electric field given by the FROG algorithm and our reconstruction from the spectrum')
% saveas(gcf,'comparison_E_fields','jpg')