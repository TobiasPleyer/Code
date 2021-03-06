function compensation_plot2(w0,t_F,Ek_F,t_opt,E_opt,t,E,w_Sk,I_Sk,p_Sk,filtered_p_Sk,P,P_opt)
    global figNum
figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(t_F,abs(Ek_F).^2,'b')
    plot(t,abs(E).^2,'g')
    % Shift it correctly
    t0 = find_closest_idx(t_opt,0);
    E0 = find_closest_idx(abs(E_opt).^2,max(abs(E_opt).^2));
    plot(t_opt,circshift(abs(E_opt).^2,t0-E0),'c')
    hold off
    title('Achieved compression compared to Fourier limit')
    xlabel('Time[fs]')
    xlim([-2000 2000])
    legend('Fourier limit','Least square','Optimization')
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
end