load('E:/Writing/Images/Matlab/Broad_Variables_Comparison_1M_Full')

plot(wavelength,P_full)
hold on
plot(wavelength,P_1M,'r-')
hold off
title('Comparison of the compressor phase compared to the phase of one mirror * N')
xlabel('wavelength [nm]','FontSize',14)
ylabel('phase [rad]','FontSize',14)
legend('Full compressor','1 mirror * N')