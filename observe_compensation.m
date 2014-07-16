%%

parent = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Ek_CHANGED/';
filebase = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
filebase2 = '26_AIR_FROG_31-W_RTT=7-us_Ip=12-A-bin';
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
             
fourier_file = sprintf('%s%s.Et-%s.dat',parent,filebase,'fourier');
Et_fourier = dlmread(fourier_file);
int_fourier = trapz(Et_fourier(:,1),Et_fourier(:,2));
trebino_file = sprintf('%s%s.Ek-%s.dat',parent,filebase,'trebino');
Ek_trebino = dlmread(trebino_file);

%%

for i=2
    fprintf('Start...\n')
    figure(1)
%%
    subplot(3,2,1)
    filename_Et = sprintf('%s%s.Et-%s.dat',parent,filebase,endings{i});
    Et = dlmread(filename_Et);
    x = Et_fourier(:,1);
    y = Et_fourier(:,2);
    plot(x,y,'color','red');
    title('Comparison between the compensated pulse and the fourier limit')
    text(500,0.75,sprintf('Curves have the\nsame integral'),'FontSize',5,'BackgroundColor',[.7 .9 .7])
    hold 'on'
    x = Et(:,1);
    y = Et(:,2);
    int = trapz(x,y);
    factor = int_fourier / int;
    y = y .* factor;
    plot(x,y,'color','green');
    hold 'off'
    axis([-2000 2000 0 1])
    %legend('Fourier limit',endings{i},'Location','SouthOutside')


    %%

    filename_Ek = sprintf('%s%s.Ek-%s.dat',parent,filebase,endings{i});
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

    subplot(3,2,2)
    [ax,h1,h2] = plotyy([Ek_trebino(:,1) Ek(:,1)],[Ek_trebino(:,2) Ek(:,2)],[Ek_trebino(:,1) Ek(:,1)],[Ek_trebino(:,3) Ek(:,3)]);
    %legend([h1;h2],'intensity trebino','intensity modified','phase trebino', 'phase modified','Location','SouthOutside')
    xlabel('time in fs')
    title('Comparison between original pulse in time and the modified one')

    gnu_file_dat_1 = sprintf('%s-%s.A',filebase,endings{i});
    fullfilename = sprintf('Daten/Dat/%s_%s_gnuplot.dat',type,gnu_file_dat_1);
    savegpbin_moritz(t, f , spec', fullfilename);

    %write gnuplot file (measured)
    gnu_file_plt=sprintf('Daten/Frog/%s_%s.plt',type,gnu_file_dat_1);
    fid_gnu_1 = fopen(gnu_file_plt,'wt');
    fprintf(fid_gnu_1,['set term jpeg\n'...
            'set output ''%sPictures/%s/%s.jpg''\n'],'Daten/',type,filebase2);
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


    %%
    name = sprintf('Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck-%s.dat',endings{i});
    Sk = dlmread(name);
    wavelength = Sk(:,1);
    phase = Sk(:,3) + Sk(:,4);
    new_phase = Sk(:,3);

    subplot(3,2,3)
    plot(wavelength,phase,'color','red')
    title('View of the whole phase curve and the smoothened area')
    xlabel('wavelength in nm')
    ylabel('phase value')
    axis([1018 1040 -15 15])
    hold 'on'
    plot(wavelength,new_phase)
    temp = (Sk(:,2).*28)-15;
    plot(wavelength,temp,'color','green')
    %legend('original','modified','Location','NorthWest','Location','SouthOutside')
    hold 'off'

    % For visual verification how the original and the modified phase differ
    subplot(3,2,5)
    plot(wavelength,phase-new_phase,'color','red')
    axis([1015 1040 -5 5])
    title('Absolute difference between original and modified phase')
    xlabel('wavelength in nm')
    ylabel('absolute difference')

    %%

    subplot(3,2,[4 6])
    file = sprintf('Daten/Pictures/%s/%s.jpg',type,filebase2);
    imshow(imread(sprintf('%s',file)));

    %%
    saveas(gcf,sprintf('Daten/Pictures/sum_up/comparison_to_fourier_%s',endings{i}),'jpg')
end