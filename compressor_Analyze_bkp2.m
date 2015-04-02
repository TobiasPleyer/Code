function compressor_Analyze()
    filename_broad_before       = 'T:\Tobias\Chirped Mirror Compressor\Analysis\004_SHG-FROG_Pure_Menlo_No_Mirrors';
    filename_broad_after_1Mirror = 'T:\Tobias\Chirped Mirror Compressor\Analysis\002_SHG-FROG_HD535_Broad_1Mirror_Bounce_Menlo';
    filename_broad_after_Full   = 'T:\Tobias\Chirped Mirror Compressor\Analysis\009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo';
    filename_broad_Reference    = 'T:\Tobias\Chirped Mirror Compressor\Analysis\001_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo_HR4000_Reference_Spectrum.txt';
    filename_narrow_before      = 'T:\Tobias\Chirped Mirror Compressor\Analysis\004_SHG_FROG_Menlo_Reference_Output (no mirrors)';
    filename_narrow_after_1Mirror = 'T:\Tobias\Chirped Mirror Compressor\Analysis\002_SHG_FROG_Compressor_Narrow_Front_Mirror_Two_Bounces';
    filename_narrow_after       = 'T:\Tobias\Chirped Mirror Compressor\Analysis\004_SHG-FROG_Narrow_Band_Compressor_Whole_New_Menlo_Settings';
    filename_narrow_Reference   = 'T:\Tobias\Chirped Mirror Compressor\Analysis\Reference_Menlo_Narrow_Compressor_Whole.txt';
    filename_narrow_sent        = 'T:\Tobias\Chirped Mirror Compressor\Originals\007_FROG_AIR_Without_multipass_RTT3.645';
    filename_broad_sent         = 'T:\Tobias\Chirped Mirror Compressor\Originals\001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A';
    filename_narrow_Trubetskov  = 'T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Narrow_Band.txt';
    filename_broad_Trubetskov   = 'T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Broad_Band.txt';
    filename_broad_Trubetskov2  = 'T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Broad_Band_no_arbitrary_constants.txt';
    
%     do_output(filename_broad_before,filename_broad_after,filename_broad_sent,filename_broad_Trubetskov2);
    %do_output(filename_broad_before,filename_broad_after_Full,filename_broad_sent,filename_broad_Trubetskov,8,8);
    do_output(filename_broad_before,filename_broad_after_1Mirror,filename_broad_sent,filename_broad_Trubetskov,8,1,true);
    %do_output(filename_narrow_before,filename_narrow_after,filename_narrow_sent,filename_narrow_Trubetskov,20,20);
    %do_output(filename_narrow_before,filename_narrow_after_1Mirror,filename_narrow_sent,filename_narrow_Trubetskov,20,2);
end

function do_output(filename_before,filename_after,filename_orig,filename_theo,N_comp,N_meas,varargin)

    if (nargin > 6)
        shift = varargin{1};
    else
        shift = false;
    end
    [pulse_info_before,...
     pulse_info_after,...
     common_Speck_before,...
     common_Speck_after,...
     spectral_phase_diff,...
     pulse_info_original,...
     Compression_Info_Measured,...
     Compression_Info_Theoretical,...
     spectral_Phase_Theoretical,...
     Int_F] = do_calculate(filename_before,filename_after,filename_orig,filename_theo,N_comp,N_meas,shift);
    
    do_plot (pulse_info_before,...
             pulse_info_after,...
             common_Speck_before,...
             common_Speck_after,...
             spectral_phase_diff,...
             pulse_info_original,...
             Compression_Info_Measured,...
             Compression_Info_Theoretical,...
             spectral_Phase_Theoretical,...
             Int_F,...
             'nosave');
    
    
%     ofile = '';
%     do_export(ofile,pulse_info);
end

function [pulse_info_before,...
          pulse_info_after,...
          common_Speck_before,...
          common_Speck_after,...
          spectral_phase_diff,...
          pulse_info_original,...
          Compression_Info_Measured,...
          Compression_Info_Theoretical,...
          spectral_Phase_Theoretical,...
          Compression_Info_Fourier] ...
    = do_calculate(filename_before,filename_after,filename_orig,filename_theo,N_comp,N_meas,shift)
    % Function do_calculate
    % This function is meant to load the stored data and bring it in the
    % correct form. Its output is meant to provide all necessary
    % information for further processing.
    postfix_Ek      = '.bin.Ek.dat';
    postfix_Speck   = '.bin.Speck.dat';
    
    file_Ek_before      = [filename_before postfix_Ek];
    file_Ek_after       = [filename_after postfix_Ek];
    file_Speck_before   = [filename_before postfix_Speck];
    file_Speck_after    = [filename_after postfix_Speck];
    
    Ek_before       = dlmread(file_Ek_before);
    Ek_after        = dlmread(file_Ek_after);
    Speck_before	= dlmread(file_Speck_before);
    Speck_after 	= dlmread(file_Speck_after);
    
    % Some retrievals have a fftshift that we have to reverse
    if shift
        Speck_after(:,2) = fftshift(Speck_after(:,2));
        Speck_after(:,3) = unwrap(fftshift(Speck_after(:,3)));
    end
    % If the FROG retrieval has a Fourier shift we have to perform those
    % additional steps in order to remove it. If this is not the case just
    % comment these lines out.
%     Speck_after(:,2) = fftshift(Speck_after(:,2));
%     Speck_after(:,3) = unwrap(fftshift(Speck_after(:,3)));
%     Ek_after(:,2) = fftshift(Ek_after(:,2));
%     Ek_after(:,3) = unwrap(fftshift(Ek_after(:,3)));
    
    pulse_info_before   = [Speck_before(:,1:3),Ek_before(:,1:3)];
    pulse_info_after    = [Speck_after(:,1:3),Ek_after(:,1:3)];
    
    pulse_info_before   = correct_pulse(pulse_info_before,'pos');
    pulse_info_after    = correct_pulse(pulse_info_after,'pos');
    
    common_wavelength = (1015:0.1:1045)';
    spectral_intensity_before   = interp1(pulse_info_before(:,1),pulse_info_before(:,2),common_wavelength);
    spectral_intensity_after    = interp1(pulse_info_after(:,1) ,pulse_info_after(:,2) ,common_wavelength);
    spectral_phase_before       = interp1(pulse_info_before(:,1),pulse_info_before(:,3),common_wavelength);
    spectral_phase_after        = interp1(pulse_info_after(:,1) ,pulse_info_after(:,3) ,common_wavelength);
    spectral_phase_diff         = spectral_phase_after-spectral_phase_before;
    spectral_phase_diff         = spectral_phase_diff ./ N_meas .* N_comp; % N_meas = # of mirrors measured
    % This is meant to accomadate several measurement cases. Example:
    % We have a compressor comprised of 20 mirrors (double pass).
    % We measure one of those mirrors but with two bounces off it.
    % Then we have N_comp = 20 and N_meas = 2
    % Phase/N_meas yields the phase of one mirror an thus
    % Phase/N_meas*N_comp gives the phase of the whole compressor
    common_Speck_before = [common_wavelength,spectral_intensity_before,spectral_phase_before];
    common_Speck_after = [common_wavelength,spectral_intensity_after,spectral_phase_after];
    [~,t_calc,Ek_calc] = compressor_toTimeDomain(spectral_intensity_before,common_wavelength,spectral_phase_before + spectral_phase_diff);
    Int_info = [t_calc,Ek_calc];
    
    % Load the original pulse and apply the measured and theoretical phase
    % to it.
    file_Ek_original     = [filename_orig postfix_Ek];
    file_Speck_original  = [filename_orig postfix_Speck];
    
    Ek_original         = dlmread(file_Ek_original);
    Speck_original      = dlmread(file_Speck_original);
    pulse_info_original = [Speck_original(:,1:3),Ek_original(:,1:3)];
    pulse_info_original = correct_pulse(pulse_info_original,'pos');
    spectral_intensity_original         = interp1(pulse_info_original(:,1),pulse_info_original(:,2),common_wavelength);
    spectral_phase_original             = interp1(pulse_info_original(:,1),pulse_info_original(:,3),common_wavelength);
    [Int_comp,t_comp,Ek_comp] = compressor_toTimeDomain(spectral_intensity_original,common_wavelength,spectral_phase_original + spectral_phase_diff,5);
    Compression_Info_Measured = {Int_comp,t_comp,Ek_comp};
    
    % Load the theoretical data
    temp = dlmread(filename_theo,'\t',1,0);
    spectral_phase_theoretical = interp1(temp(:,1),temp(:,2),common_wavelength);
    spectral_GD_theoretical = interp1(temp(:,1),temp(:,3),common_wavelength);
    spectral_GDD_theoretical = interp1(temp(:,1),temp(:,4),common_wavelength);
    clear temp
    spectral_phase_theoretical(isnan(spectral_phase_theoretical)) = 0;
    spectral_GD_theoretical(isnan(spectral_GD_theoretical)) = 0;
    spectral_GDD_theoretical(isnan(spectral_GDD_theoretical)) = 0;
    spectral_phase_theoretical = N_comp .* spectral_phase_theoretical; % This is necessary because the data was provided as per mirror, N_comp bounces in the compressor
    spectral_GD_theoretical = N_comp.* spectral_GD_theoretical; % N_comp = # of mirrors in the compressor
    spectral_GDD_theoretical = N_comp.* spectral_GDD_theoretical;
    spectral_Phase_Theoretical = {spectral_phase_theoretical,spectral_GD_theoretical,spectral_GDD_theoretical};
    [Int_theo,t_theo,Ek_theo] = compressor_toTimeDomain(spectral_intensity_original,common_wavelength,spectral_phase_original - spectral_phase_theoretical,5);
    Compression_Info_Theoretical1 = {Int_theo,t_theo,Ek_theo};
    [Int_theo,t_theo,Ek_theo] = compressor_toTimeDomain(spectral_intensity_original,common_wavelength,spectral_phase_before + spectral_phase_theoretical,5);
    Compression_Info_Theoretical2 = {Int_theo,t_theo,Ek_theo};
    [Int_theo,t_theo,Ek_theo] = compressor_toTimeDomain(spectral_intensity_original,common_wavelength,spectral_phase_before - spectral_phase_theoretical,5);
    Compression_Info_Theoretical3 = {Int_theo,t_theo,Ek_theo};
    Compression_Info_Theoretical = {Compression_Info_Theoretical1,Compression_Info_Theoretical2,Compression_Info_Theoretical3};
    
    % Provide the Fourier limit of the spectrum in order to scale the
    % graphs correctly
    [Int_F,t_F,Ek_F] = compressor_toTimeDomain(spectral_intensity_original,common_wavelength,3);
    Compression_Info_Fourier = {Int_F,t_F,Ek_F};
end

function do_plot(pulse_info_before,...
                 pulse_info_after,...
                 common_Speck_before,...
                 common_Speck_after,...
                 Phase_Diff,...
                 pulse_info_original,...
                 Compression_Info_Measured,...
                 Compression_Info_Theoretical,...
                 spectral_Phase_Theoretical,...
                 Compression_Info_Fourier,...
                 fsave)
    % Function do_plot
    % This function generates all plots to viauslize the data.
    % The arguments are vectors and matrices that include all necessary
    % data in column form.
    %       e.g.: retr_Before = [l,I_l,p_l,t,I_t,p_t]
    styles = {'r-', 'b-', 'g-', 'c-', 'm-', 'k-'};
    colors = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};
    F_size_Label = 22;
    % Plot the retrieved spectra before and after the compressor plus the
    % reference spectrum and the resulting phase of the compressor
    full_wavelength_before   = pulse_info_before(:,1);
    wintensity_before   = pulse_info_before(:,2);
    wphase_before       = pulse_info_before(:,3);
    time_before         = pulse_info_before(:,4);
    tintensity_before   = pulse_info_before(:,5);
    tphase_before       = pulse_info_before(:,6);
    full_wavelength_after    = pulse_info_after(:,1);
    wintensity_after    = pulse_info_after(:,2);
    wphase_after        = pulse_info_after(:,3);
    time_after          = pulse_info_after(:,4);
    tintensity_after    = pulse_info_after(:,5);
    tphase_after        = pulse_info_after(:,6);
    wavelength_before   = common_Speck_before(:,1);
    intensity_before    = common_Speck_before(:,2);
    phase_before        = common_Speck_before(:,3);
    wavelength_after    = common_Speck_after(:,1);
    intensity_after     = common_Speck_after(:,2);
    phase_after         = common_Speck_after(:,3);
    Int_theo2       = Compression_Info_Theoretical{2}{1};
    time_theo2      = Compression_Info_Theoretical{2}{2};
    Ek_theo2        = Compression_Info_Theoretical{2}{3};
    Int_theo3       = Compression_Info_Theoretical{3}{1};
    time_theo3      = Compression_Info_Theoretical{3}{2};
    Ek_theo3        = Compression_Info_Theoretical{3}{3};
    spectral_phase_theoretical = spectral_Phase_Theoretical{1};
    phase_expected  = phase_before + spectral_phase_theoretical;
    %phase_expected = tilt_curve(wavelength_before,phase_after,phase_expected);
    figure()
        subplot(2,2,1)
            H2 = zeros(1,4);
            [AX, H1, H2(1)] = plotyy(full_wavelength_before,wintensity_before,full_wavelength_before,wphase_before);
            set(H1,'Color',colors{1},'LineStyle','-')
            set(H2(1),'Color',colors{1},'LineStyle','-')
            hold(AX(2))
            [H2(2)] = plot(AX(2),wavelength_before,spectral_phase_theoretical);
            set(H2(2),'Color',colors{2},'LineStyle','--')
            [H2(3)] = plot(AX(2),wavelength_before,phase_after + phase_before);
            set(H2(3),'Color',colors{3},'LineStyle','--')
            [H2(4)] = plot(AX(2),wavelength_before,phase_after - phase_before);
            set(H2(4),'Color',colors{4},'LineStyle','--')
            title(sprintf('Menlo FROG spectrum before compressor'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            xlim(AX(1),[1000 1060])
            xlim(AX(2),[1000 1060])
            ylim(AX(2),[-100 150])
            l = legend(H2,{'P_{Measured}', 'P_{Theoretical}', 'P_{after} + P_{before}', 'P_{after} - P_{before}'},'Location','NorthWest');
        subplot(2,2,2)
            H2 = zeros(1,3);
            [AX,H1,H2(1)] = plotyy(wavelength_after,intensity_after,wavelength_after,phase_after);
            set(H1,'Color',colors{1},'LineStyle','-')
            set(H2(1),'Color',colors{1},'LineStyle','-')
            hold(AX(1))
            hold(AX(2))
            [H2(2)] = plot(AX(2),wavelength_before,phase_expected);
            set(H2(2),'Color',colors{4},'LineStyle','--')
            [H2(3)] = plot(AX(2),wavelength_before,phase_before - spectral_phase_theoretical);
            set(H2(3),'Color',colors{5},'LineStyle','--')
            title(sprintf('Menlo FROG spectrum after compressor'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            %xlim(AX(1),[1014 1045])
            %xlim(AX(2),[1014 1045])
            ylim(AX(2),[-200 100])
            l = legend(H2,{'P_{Measured}', 'P_{before} + P_{theo}', 'P_{before} - P_{theo}'},'Location','NorthWest');
            %set(l,'FontSize',F_size_Label);
        subplot(2,2,3)
            [AX, H1, H2] = plotyy(time_before,tintensity_before,time_before,tphase_before);
            title(sprintf('Menlo FROG temporal measurement before compressor'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
        subplot(2,2,4)
            H1 = zeros(1,3);
            H2 = zeros(1,3);
            [AX, H1(1), H2(1)] = plotyy(time_after,tintensity_after,time_after,tphase_after);
            set(H1(1),'Color',colors{1},'LineStyle','-')
            set(H2(1),'Color',colors{1},'LineStyle','--')
            hold(AX(1))
            hold(AX(2))
            [H1(2)] = plot(AX(1),time_theo2,abs(Ek_theo2).^2);
            [H2(2)] = plot(AX(2),time_theo2,unwrap(-angle(Ek_theo2)));
            set(H1(2),'Color',colors{2},'LineStyle','-')
            set(H2(2),'Color',colors{2},'LineStyle','--')
            [H1(3)] = plot(AX(1),time_theo3,abs(Ek_theo3).^2);
            [H2(3)] = plot(AX(2),time_theo3,unwrap(-angle(Ek_theo3)));
            set(H1(3),'Color',colors{3},'LineStyle','-')
            set(H2(3),'Color',colors{3},'LineStyle','--')
            title(sprintf('Menlo FROG temporal measurement after compressor'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            l = legend(H1,{'Measured', 'Meas + Theo', 'Meas - Theo'},'Location','NorthWest');
%         subplot(2,2,4)
%             H1 = zeros(1,2); % Allocate the array for the handles
%             H2 = zeros(1,2);
%             [AX,H1(1),H2(1)] = plotyy();
%             title('Fourier transformations of several phase combinations','FontSize',14)
%             set(AX(1),'YColor','black')
%             set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
%             set(get(AX(1),'YLabel'),'String','Intensity/Integral','FontSize',F_size_Label)
%             set(AX(2),'YColor','black')
%             set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
%             set(H1(1),'Color',colors{1},'LineWidth',2)
%             set(H2(1),'Color',colors{1},'LineStyle','--')
%             hold(AX(1))
%             hold(AX(2))
%             [H1(2)] = plot(AX(1),t_2,abs(Ek_2).^2./Int_2);
%             [H2(2)] = plot(AX(2),t_2,unwrap(-angle(Ek_2)));
%             set(H1(2),'Color',colors{2},'LineWidth',2)
%             set(H2(2),'Color',colors{2},'LineStyle','--')
%             %ylim(AX(1),[0 0.0013])
%             %xlim(AX(1),[-10000 10000])
%             %xlim(AX(2),[-10000 10000])
%             l = legend(H1,{'V1', 'V2', 'V3'},'Location','NorthWest');
%             set(l,'FontSize',F_size_Label);
    % Plot the calculated GD and GDD values for the measured phase
    [P,GD,GDD] = GDD_from(wavelength_before,Phase_Diff);
    GD_theoretical = spectral_Phase_Theoretical{2};
    GDD_theoretical = spectral_Phase_Theoretical{3};
    %P_theo = tilt_curve(wavelength_before,Phase_Diff,spectral_phase_theoretical);
    figure()
        subplot(2,2,1)
            leg = {};
            plot(wavelength_before,Phase_Diff,'r-','LineWidth',2)
            leg{1} = sprintf('Measured Phase');
            hold on
            plot(wavelength_before,P,'g-','LineWidth',2)
            leg{2} = sprintf('Polynomial Fit');
            plot(wavelength_before,spectral_phase_theoretical,'c-','LineWidth',2)
            leg{3} = sprintf('Theoretical Phase');
            hold off
            xlim([1000 1060])
            title(sprintf('Measured and fitted Phase Function\nof the Chirped Mirror Compressor'),'FontSize',22)
            xlabel('wavelength [nm]','FontSize',F_size_Label)
            ylabel('phase [rad]','FontSize',F_size_Label)
            l = legend(leg,'Location','NorthWest');
            set(l,'FontSize',F_size_Label);
        subplot(2,2,2)
            % Before we can plot the data we have to make sure that they have the same
            % offset from the axis in order to be comparable. The value of the
            % theoretical measurement is taken as zero point
            plot(wavelength_before,GD,'g-')
            hold on
            plot(wavelength_before,GD_theoretical,'c-','LineWidth',2)
            hold off
            xlim([1000 1060])
            title(sprintf('Measured GD for the\nchirped mirror compressor'),'FontSize',22)
            xlabel('wavelength [nm]','FontSize',F_size_Label)
            ylabel('GD [fs]','FontSize',F_size_Label)
            l = legend('Measured GD','Theoretical GD','Location','SouthEast');
            set(l,'FontSize',F_size_Label);
        subplot(2,2,3)
            plot(wavelength_before,GDD,'g-')
            hold on
            plot(wavelength_before,GDD_theoretical,'c-','LineWidth',2)
            hold off
            xlim([1000 1060])
            title(sprintf('Measured GDD for the\nchirped mirror compressor'),'FontSize',22)
            xlabel('wavelength [nm]','FontSize',F_size_Label)
            ylabel('GDD [fs^2]','FontSize',F_size_Label)
            l = legend('Measured GDD','Theoretical GDD','Location','SouthEast');
            set(l,'FontSize',F_size_Label);
    
    % Plot the original pulse in time domain. Plot what we would expect
    % from the perfect theoretical phase. Plot the retrieved pulse in time
    % domain.
    wavelength_orig = pulse_info_original(:,1);
    wintensity_orig = pulse_info_original(:,2);
    wphase_orig     = pulse_info_original(:,3);
    time_orig       = pulse_info_original(:,4);
    tintensity_orig = pulse_info_original(:,5);
    tphase_orig     = pulse_info_original(:,6);
    Int_comp        = Compression_Info_Measured{1};
    time_comp       = Compression_Info_Measured{2};
    Ek_comp         = Compression_Info_Measured{3};
    Int_theo1       = Compression_Info_Theoretical{1}{1};
    time_theo1      = Compression_Info_Theoretical{1}{2};
    Ek_theo1        = Compression_Info_Theoretical{1}{3};
    Int_Orig        = abs(trapz(time_orig,tintensity_orig));
    Int_F           = Compression_Info_Fourier{1};
    t_F             = Compression_Info_Fourier{2};
    Ek_F            = Compression_Info_Fourier{3};
    figure()
        subplot(2,2,1)
            [AX, H1, H2] = plotyy(time_orig,tintensity_orig,time_orig,tphase_orig);
            title(sprintf('The original pulse provided for the chirped mirror production'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            ylim(AX(1),[0 1])
            xlim(AX(1),[-0.4e4 0.4e4])
            xlim(AX(2),[-0.4e4 0.4e4])
        subplot(2,2,2)
            H1 = zeros(1,2);
            [AX, H1(1), H2] = plotyy(time_theo1,abs(Ek_theo1).^2./Int_theo1.*Int_F,time_theo1,unwrap(-angle(Ek_theo1)));
            set(H1(1),'Color',colors{1},'LineStyle','-')
            set(H2,'Color',colors{3},'LineStyle','-')
            hold(AX(1))
            [H1(2)] = plot(AX(1),t_F,abs(Ek_F).^2);
            set(H1(2),'Color',colors{2},'LineStyle','--')
            title(sprintf('The compression result we expect from our theoretical phase'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            ylim(AX(1),[0 1])
            xlim(AX(1),[-0.2e4 0.3e4])
            xlim(AX(2),[-0.2e4 0.3e4])
        subplot(2,2,3)
            H1 = zeros(1,2);
            [AX, H1(1), H2] = plotyy(time_comp,abs(Ek_comp).^2./Int_comp.*Int_F,time_comp,unwrap(-angle(Ek_comp)));
            set(H1(1),'Color',colors{1},'LineStyle','-')
            set(H2,'Color',colors{3},'LineStyle','--')
            hold(AX(1))
            [H1(2)] = plot(AX(1),t_F,abs(Ek_F).^2);
            set(H1(2),'Color',colors{2},'LineStyle','--')
            title(sprintf('The compression result we expect from our measured phase'),'FontSize',14)
            set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
            set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
            set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
            ylim(AX(1),[0 1])
            xlim(AX(1),[-1e4 1e4])
            xlim(AX(2),[-1e4 1e4])

%     quick_Test_for_GDD_Calc()
    fprintf('Successfully plotted.\n');
    if (strcmp(fsave,'save'))
        %...
        fprintf('Successfully saved plots.\n')
    end
end

function do_export(ofile,pulse_info)
    % Function do_export
    % Writes the given pulse information in two seperate files.
    % The format and naming conventions of these files are very similar to
    % those of Rick Trebino and his team.
    % The argument is expected to be given in the following form:
    %       pulse_info = [lambda,I_lambda,p_lambda,time,I_time,p_time]
    Ek  = [ofile '.Ek.dat'];
    Sk  = [ofile '.Speck.dat'];
    l   = pulse_info(:,1);
    I_l = pulse_info(:,2);
    p_l = pulse_info(:,3);
    t   = pulse_info(:,4);
    I_t = pulse_info(:,5);
    p_t = pulse_info(:,6);
    %Prepare spectral pulse information for export
    M = [l,I_l,p_l];
    dlmwrite(Ek,M,'\t');
    %Prepare temporal pulse information for export
    M = [t,I_t,p_t];
    dlmwrite(Sk,M,'\t');
    fprintf('Files exported.\n')
end

function cpulse_info=correct_pulse(pulse_info,varargin)
    % Funciton correct_phase
    % This function makes sure the GDD of the phase has the desired sign as
    % indicated by the second argument. The default is positive GDD.
    % p = correct_phase(pulse_info[,'pos'/'neg']);
    if (nargin == 1)
        GDD = 'pos';
    else
        GDD = varargin{1};
    end
    phase = pulse_info(:,3); %this is supposed to be p_t
    L = length(phase);
    [~,I] = min(phase);
    if (I(1) < 20 || I(1) > (L-20)) %Not the most sophisticated check. Only works when the GDD is more dominant then the higher order, which is the case for most of our cases
        sign = 'neg';
    else
        sign = 'pos';
    end
    if (strcmp(GDD,sign))
        cpulse_info = pulse_info;
    else
        l   = pulse_info(:,1);
        I_l = pulse_info(:,2);
        p_l = pulse_info(:,3);
        t   = pulse_info(:,4);
        I_t = pulse_info(:,5);
        p_t = pulse_info(:,6);
        p_l = -p_l;
        I_t = fliplr(I_t')';
        p_t = -fliplr(p_t')';
        cpulse_info = [l,I_l,p_l,t,I_t,p_t];
    end
end

function [P,GD,GDD]=GDD_from(wavelength,phase)
    % Function GDD_from
    % This function aims to calculate the GD and GDD values of the provided
    % phase. This is done in the following steps:
    %       1. Transform into frequency domain
    %       2. Fit a polynomial to the phase
    %       3. In frequency domain form the derivatives of that polynomial
    %       4. GD and GDD are the 1st and 2nd of those derivatives
    c = 299792458;
    % We know the phase in wavelength domain -> transform to frequency space
    % and calculate the GD and GDD
%     frequency = (2*pi*c) ./ (wavelength .* 1e-9);
%     poly = polyfit(frequency,-phase,2);
%     P = polyval(poly,frequency);
%     poly_deriv = polyder(poly);
%     GD = polyval(poly_deriv,frequency) .* 1e15;
%     poly_deriv = polyder(poly_deriv);
%     GDD = polyval(poly_deriv,frequency) .* 1e30;
    frequency = (2*pi*c) ./ (wavelength .* 1e-9);
    frequency = fliplr(frequency')'; % The fliplr comes because we change our x-axis: l --> w => ascending --> descending
    phase = fliplr(phase')';
    poly = polyfit(frequency,phase,6);
    P = polyval(poly,frequency);
    poly_deriv = polyder(poly);
    GD = polyval(poly_deriv,frequency) .* 1e15;
    poly_deriv = polyder(poly_deriv);
    GDD = polyval(poly_deriv,frequency) .* 1e30;
    P = fliplr(P')';
    GD = fliplr(GD')';
    GDD = fliplr(GDD')';
end

function prove_validity()
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
    file_name = 'FROG/009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo.bin.Speck.dat';
    S = dlmread(file_name);
    wavelength_with_mirrors = S(:,1);
    spectral_intensity_with_mirrors = S(:,2);
    spectral_phase_with_mirrors = S(:,3);
    file_name = 'FROG/009_SHG-FROG_HD535_Broad_Chirped_Mirror_Compressor_Whole_Menlo.bin.Ek.dat';
    S = dlmread(file_name);
    time_with_mirrors = S(:,1);
    temporal_intensity_with_mirrors = S(:,2);
    temporal_phase_with_mirrors = S(:,3);
    figure()
    subplot(2,2,1)
    plotyy(wavelength_with_mirrors,spectral_intensity_with_mirrors,wavelength_with_mirrors,spectral_phase_with_mirrors)
    subplot(2,2,2)
    plotyy(time_with_mirrors,temporal_intensity_with_mirrors,time_with_mirrors,temporal_phase_with_mirrors)
    subplot(2,2,3)
    [Int,t,E] = compressor_toTimeDomain(spectral_intensity_with_mirrors,wavelength_with_mirrors,spectral_phase_with_mirrors);
    plotyy(t,abs(E).^2,t,unwrap(-angle(E)))
    subplot(2,2,4)
    plotyy(t,abs(E).^2-temporal_intensity_with_mirrors',t,unwrap(-angle(E))-temporal_phase_with_mirrors')
end

function quick_Test_for_GDD_Calc()
    B = dlmread('T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Broad_Band.txt','\t',1,0);
    l_T = B(:,1);
    P_T = B(:,2);
    GD_T = B(:,3);
    GDD_T = B(:,4);
    
    [P,GD,GDD] = GDD_from(l_T,P_T);
    
    figure()
    subplot(1,2,1)
    plot(l_T,-GD)
    hold on
    plot(l_T,GD_T,'r-')
    hold off
    xlabel('wavelength [nm]')
    ylabel('GD [fs]')
    legend('Pleyer','Trubetskov','Location','NorthWest')
    subplot(1,2,2)
    plot(l_T,-GDD)
    hold on
    plot(l_T,GDD_T,'r-')
    hold off
    xlabel('wavelength [nm]')
    ylabel('GDD [fs^2]')
    legend('Pleyer','Trubetskov','Location','NorthWest')
end

function [p2_out]=tilt_curve(w,p1_in,p2_in)
    % Function tile_curve
    % This function aims to add a linear term to an existing phase curve in
    % order to be better comparable to another one.
    % This is no optimization method but simply uses a mere geometrical
    % idea. The returned phase p2_out will have the same GDD and higher
    % dispersion coefficients as p2_in, but with linear and constant term
    % closer to those of p1_in.
    % It is expected that w is a wavelength range common to both input
    % phase curves!
    D = p2_in-p1_in;
    idx = auxiliary_find_closest_idx(D,0);
    % Make sure we are not too far to the borders of the array
    if (idx <= 100 || length(w)-idx >= 100)
        m = min(idx,length(w)-idx);
        idx1 = idx - m + 2;
        idx2 = idx + m - 2;
    else
        idx1 = idx - 100;
        idx2 = idx + 100;
    end
    x1 = w(idx1);
    x2 = w(idx2);
    y1 = p2_in(idx1) - p1_in(idx1);
    y2 = p2_in(idx2) - p1_in(idx2);
    %p2_out = p2_in + ((y2 - y1)/(x2-x1)).*w + (y1*x2 - y2*x1)/(x2-x1); 
    p2_out = p2_in - ((y2 - y1)/(x2-x1)).*(w-w(idx));
end