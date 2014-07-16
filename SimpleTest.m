clear;

dir_main = 'Daten\';
filename_base = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';

filename_Ek = sprintf('%s%s.Ek.dat',dir_main, filename_base);
Ek_orig = dlmread(filename_Ek);
Ek_orig = Ek_orig(:,1:3);
% filename_rec = sprintf('%s%s.Arecon.dat',dir_main, filename_base);
% recon = dlmread(filename_rec);

lambda0 = 1030;
l_min = 1010;
l_max = 1045;
gridSize = 1000;
tau_range = 5000;
tbuf = 000;

[t, f, spec, type] = PG_FROG(Ek_orig, tau_range, tbuf, lambda0, l_min, l_max, gridSize);

[TAU, F] = meshgrid(t, f);


figure(1)
surf(TAU, F, spec,'EdgeColor','none');
view(0,90)
title(type)
xlabel('Time [fs]')
ylabel('Frequency [1/fs]')

    gnu_file_dat_1 = sprintf('%s.A',filename_base);
    fullfilename = sprintf('%s%s_%s_gnuplot.dat',dir_main,type,gnu_file_dat_1);
    savegpbin_moritz(t, f , spec', fullfilename);

    %write gnuplot file (measured)
    gnu_file_plt=sprintf('%s%s_%s.plt',dir_main,type,gnu_file_dat_1);
    fid_gnu_1 = fopen(gnu_file_plt,'wt');
    fprintf(fid_gnu_1,['set term jpeg\n'...
        'set output ''%s%s_%s.jpeg''\n'],dir_main,type,gnu_file_dat_1);
    str = sprintf(['set pm3d map\n'...
        'set palette defined (0 "black", 0.01 "red",0.16 "yellow",0.46 "green", 0.81 "blue", 1 "white")\n'...
        'set title ''%s trace\n'...
        'set ytics auto\n'...
        'set cblabel ''Intensity (a.u)''\n'...
        'set xtics auto\n'...
        'set yrange [:]\n'...
        'set xlabel ''Delay [fs]''\n'...
        'set ylabel ''Wavelength [nm]''\n'...
        'set xrange [:]\n'...
        'splot ''%s'' binary matrix with pm3d notitle\n'],...
        type,fullfilename);
    fprintf(fid_gnu_1,'%s',str);
    fprintf(fid_gnu_1,'set output ');

    fclose(fid_gnu_1);

    %execute gnuplot to create jpeg
    gnu_file_plt_eval = sprintf(' "%s"',gnu_file_plt);
    eval(['!"C:\Program Files (x86)\gnuplot\bin\gnuplot.exe"', gnu_file_plt_eval]);