load('T:/Tobias/Chirped Mirror Compressor/Analysis/Broad_Variables_Comparison_1M_Full')

plot(wavelength,P_full)
hold on
plot(wavelength,P_1M,'r-')
hold off
%title('Comparison of the compressor phase compared to the phase of one mirror * N')
xlabel('wavelength [nm]','FontSize',16)
ylabel('phase [rad]','FontSize',16)
h = legend(sprintf('Full Compressor\n(8 Bounce HD535)'),'HD535 * 8');
set(gca,'FontSize',16)
set(h,'FontSize',14)
saveas(gca,'../Bilder/Comparison_1M_Full','epsc')