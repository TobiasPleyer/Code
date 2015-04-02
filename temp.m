A = dlmread('T:\Tobias\Chirped Mirror Compressor\001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A.bin.Speck.dat');
B = dlmread('T:\Tobias\Chirped Mirror Compressor\Trubetskov_Mirror_Specification_Broad_Band_All.txt','\t',1,0);
w = B(:,3);
P = B(:,4);
P = P(w>0);
w = w(w>0);

figure()
subplot(1,2,1)
plot(A(:,1),A(:,3))
hold on
plot(w,P,'r-')
hold off
xlim([1015,1045])
xlabel('wavelength [nm]')
ylabel('Phase [rad]')
legend('T. Pleyer','M. Trubetskov')
subplot(1,2,2)
plot(A(:,1),A(:,3))
hold on
plot(w,8*P,'r-')
hold off;
xlim([1015,1045])
xlabel('wavelength [nm]')
ylabel('Phase [rad]')
legend('T. Pleyer','8* M. Trubetskov')

C = dlmread('T:\Tobias\Chirped Mirror Compressor\broadband_pulse_intensity.dat',',');
D = dlmread('T:\Tobias\Chirped Mirror Compressor\broadband_pulse_phase.dat',',');

figure()
plotyy(C(:,1),C(:,2),D(:,1),D(:,2))