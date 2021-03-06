function compensation_plot(w0,t_Et,I_Et,p_Et,w_Sk,I_Sk,p_Sk,filtered_p_Sk,D1,D2,D3,D2_opt,D3_opt,P,P_opt,t_F,Ek_F,t,E,t_opt,E_opt,max_order,opt_peaks,peaks)
    global figNum
    %-------------------------Figure(1)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    plotyy(t_Et,I_Et,t_Et,p_Et)
    title('Intensity and phase of our pulse in the time domain')
    xlabel('Time[fs]')
    legend('Intensity','Phase')
    %-------------------------Figure(2)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    plotyy(w_Sk,I_Sk,w_Sk,p_Sk)
    title('Intensity and phase of our pulse in the frequency domain')
    xlabel('Omega[1/fs]')
    legend('Intensity','Phase')
    %-------------------------Figure(3)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,p_Sk)
    plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
    hold off
    title('Phase in spectral domain')
    xlabel('Omega[1/fs]')
    %-------------------------Figure(4)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk(1:end-1),D1)
    plot([w0 w0],[min(D1) max(D1)],'r')
    hold off
    title('GD of our pulse (filter applied)')
    xlabel('Omega[1/fs]')
    %-------------------------Figure(5)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk(1:end-2),D2)
    plot(w_Sk(1:end-2),D2_opt,'g')
    plot([w0 w0],[min(D2) max(D2)],'r')
    hold off
    title('GDD of our pulse (filter applied)')
    xlabel('Omega[1/fs]')
    legend('filtered phase','custom')
    %-------------------------Figure(6)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk(1:end-3),D3)
    plot(w_Sk(1:end-3),D3_opt,'g')
    plot([w0 w0],[min(D3) max(D3)],'r')
    hold off
    title('TOD of our pulse (filter applied)')
    xlabel('Omega[1/fs]')
    legend('filtered phase','custom')
    %-------------------------Figure(7)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,I_Sk./max(I_Sk).*max(filtered_p_Sk),'k')
    plot(w_Sk,p_Sk,'b')
    plot(w_Sk,P,'r')
    plot(w_Sk,P_opt,'g')
    plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
    hold off
%     xlim([lower-0.01 higher+0.01])
    ylim([min(p_Sk)-2 max(p_Sk)+2])
    title('Taylor approximations for our phase curve')
    xlabel('Omega[1/fs]')
    legend('Intensity','Original Phase','least square','custom')
    %-------------------------Figure(8)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,I_Sk./max(I_Sk).*50,'k')
    plot(w_Sk,p_Sk-P,'g')
    plot(w_Sk,p_Sk-P_opt,'c')
    plot([w0 w0],[min(p_Sk) max(p_Sk)],'r')
    hold off
%     xlim([lower-0.01 higher+0.01])
    ylim([-20 50])
    title('Difference between the Taylor approximation and real phase')
    xlabel('Omega[1/fs]')
    legend('Intensity','Original Phase','least square','custom')
    %-------------------------Figure(9)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(t_F,abs(Ek_F).^2,'b')
    plot(t,abs(E).^2,'g')
    % Shift it correctly
    t0 = find_closest_idx(t_opt,0);
    E0 = find_closest_idx(abs(E_opt).^2,max(abs(E_opt).^2));
    plot(t_opt,circshift(abs(E_opt).^2,[1,t0-E0]),'c')
    hold off
    title('Achieved compression compared to Fourier limit')
    xlabel('Time[fs]')
    xlim([-2000 2000])
    legend('Fourier limit','Least square','Optimization')
    %-------------------------Figure(10)------------------------------------
    %----------------------------------------------------------------------
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(1:max_order,opt_peaks)
    plot(1:max_order,peaks,'r')
    hold off
    title('Compression quality with increasing order')
    xlabel('Order')
    ylabel('Peak percentage of Fourier limit')
    legend('custom','least squares')
end