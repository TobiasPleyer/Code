function [common_wavelength,Speck_int,Speck_phase]=auxiliary_plotPulse(filename,FWHM_flag,title_str,common_wavelength)
    %Function auxiliary_plotPulse(filename[FWHM,title_str,common_wavelength]).
    %Arguments:
    %           filename:          string Mandatory
    %                              This string gives the location and
    %                              filebase of the data to be plotted.
    %           FWHM:              boolean optional
    %                              If true the FWHM will be drawn into the
    %                              graph.
    %           common_wavelength: vector optional
    %                              If provided this vector will be made the
    %                              new X-axis data for the plots.
    %           title_str:         string optional
    %                              The title to be displayed on the graphs.
    %This function plots the temporal and spectral characterisitc of the
    %pulse provided by the files with the common base "filename".
    %The filename argument should be given without suffixes, these will be
    %added automatically.
    %If the optional argument "common_wavelength" is given, the function
    %will interpolate and resample the data to have "common_wavelength" as
    %the new X-axis data.
    
    Speck_title = sprintf('Spectral pulse shape of %s',filename);
    Ek_title = sprintf('Temporal pulse shape of %s',filename);
    cw = false;
    FWHM = false;
    
    if nargin == 2
        FWHM = FWHM_flag;
        cw = false;
    elseif nargin == 3
        Speck_title = sprintf('Spectral_domain_%s',title_str);
        Ek_title = sprintf('Temporal_domain_%s',title_str);
        FWHM = FWHM_flag;
    elseif nargin == 4
        Speck_title = sprintf('Spectral_domain_%s',title_str);
        Ek_title = sprintf('Temporal_domain_%s',title_str);
        FWHM = FWHM_flag;
        cw = true;
    else % Should never be entered
        error('Wrong number of arguments! See "help auxiliary_plotPulse"')
    end
    F_size_Label = 16;

    filename_cmp_Speck = [filename '.bin.Speck.dat'];
    temp = dlmread(filename_cmp_Speck);
        Speck_wavel = temp(:,1);
        Speck_int = temp(:,2);
        Speck_phase = temp(:,3);

    if cw
        Speck_int = interp1(Speck_wavel,Speck_int,common_wavelength);
        Speck_phase = interp1(Speck_wavel,Speck_phase,common_wavelength);
    else
        common_wavelength = Speck_wavel;
    end

    [~,t,Ek] = compressor_toTimeDomain(Speck_int,common_wavelength,Speck_phase,5);
    [~,t2,Ek2] = compressor_toTimeDomain(Speck_int,common_wavelength,Speck_phase);
    figure()
    [AX,~,~] = plotyy(t,abs(Ek).^2,t2,unwrap(-angle(Ek2)));
    if FWHM
        [t_min,t_max,~] = auxiliary_findFWHM(t,abs(Ek).^2);
        hold(AX(1))
        plot(AX(1),[t_min,t_min,t_max,t_max],[0,0.5,0.5,0],'r-')
        label = sprintf('FWHM: %2.0f fs',t_max-t_min);
        text((t_max+t_min)/2-length(label)*30,0.45,label)
    end
    %title(Ek_title,'FontSize',F_size_Label)
    set(get(AX(1),'XLabel'),'String','time [fs]','FontSize',F_size_Label)
    set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
    set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
    xlim(AX(1),[-3000 3000])
    xlim(AX(2),[-3000 3000])
    %The y-scale is calculated from all availbale x-data, not just the one
    %used by xlim -> find the right scaling
    idx1 = auxiliary_find_closest_idx(t2,-1500);
    idx2 = auxiliary_find_closest_idx(t2,1500);
    p =unwrap(-angle(Ek2));
    m = floor(min(p(min(idx1,idx2):max(idx1,idx2))));
    M = ceil(max(p(min(idx1,idx2):max(idx1,idx2))));
    ylim(AX(2),[m M])
    if M-m > 30
        set(AX(2),'YTick',m:4:M)
    else
        set(AX(2),'YTick',m:2:M)
    end
    set(AX(1),'Box','off')
    set(AX,'FontSize',F_size_Label)
    saveas(gcf,sprintf('../Bilder/%s',Ek_title),'epsc')

    figure()
    [AX,~,~] = plotyy(common_wavelength,Speck_int,common_wavelength,Speck_phase);
    %title(Speck_title,'FontSize',F_size_Label)
    set(get(AX(1),'XLabel'),'String','wavelength [nm]','FontSize',F_size_Label)
    set(get(AX(1),'YLabel'),'String','intensity [a.u.]','FontSize',F_size_Label)
    set(get(AX(2),'YLabel'),'String','phase [rad]','FontSize',F_size_Label)
    x = 1020;
    X = 1040;
    xlim(AX(1),[x X])
    xlim(AX(2),[x X])
    %The y-scale is calculated from all availbale x-data, not just the one
    %used by xlim -> find the right scaling
    idx1 = auxiliary_find_closest_idx(common_wavelength,x);
    idx2 = auxiliary_find_closest_idx(common_wavelength,X);
    m = floor(min(Speck_phase(min(idx1,idx2):max(idx1,idx2))));
    M = ceil(max(Speck_phase(min(idx1,idx2):max(idx1,idx2))));
    ylim(AX(2),[m M])
    if M-m > 30
        set(AX(2),'YTick',m:4:M)
    else
        set(AX(2),'YTick',m:2:M)
    end
    set(AX(1),'Box','off')
    set(AX,'FontSize',F_size_Label)
    saveas(gcf,sprintf('../Bilder/%s',Speck_title),'epsc')
end