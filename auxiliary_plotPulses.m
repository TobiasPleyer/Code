function auxiliary_plotPulses(folder_base,names,expl,titl,Int_F)
    % Function auxiliary_plotPulses(folder_base,names,expl,title[,Int_F])
    
    if nargin == 5
        F = true;
    else
        F = false;
    end
    F_size_Label = 16;
    colors = {[0 0 1], [1 0 0], [0 1 0], [0 1 1], [1 0 1], [0 0 0], [1 0.6 0]};
    
    
    figure()
    hold on
    leg = {};
    for i=1:length(names)
        filename_Speck = [folder_base names{i} '.bin.Speck.dat'];
        temp = dlmread(filename_Speck);
            Speck_wavel = temp(:,1);
            Speck_int = temp(:,2);
            Speck_phase = temp(:,3);
        [Int,t,Ek] = compressor_toTimeDomain(Speck_int,Speck_wavel,Speck_phase,5);
        if ~F
            [Int_F,~,~] = compressor_toTimeDomain(Speck_int,Speck_wavel,5);
        end
        plot(t,abs(Ek).^2./Int.*Int_F,'Color',colors{i})
        leg{i} = expl{i};
    end
    hold off
    legend(leg)
    title(titl,'FontSize',14)
    xlabel('time [fs]','FontSize',F_size_Label)
    ylabel('intensity a.u.','FontSize',F_size_Label)
    xlim([-2000 2000])
    set(gca,'FontSize',F_size_Label)
    
end