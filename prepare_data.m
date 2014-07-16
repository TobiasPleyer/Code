clear

F=dlmread('T:/Tobias/Daten/20140701_Chirped_Mirror_FROGS/001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A.bin.Speck.dat');

wavelength = min(F(:,1)):0.1:max(F(:,1));
vq=interp1(F(:,1),F(:,2),wavelength);
intensity = vq;
vq=interp1(F(:,1),F(:,3),wavelength);
phase = vq;
subplot(2,2,1)
plot(F(:,1),F(:,2))
title('Original Intensity')
subplot(2,2,2)
plot(wavelength,intensity)
title('Interpolated Intensity')
subplot(2,2,3)
plot(F(:,1),F(:,3))
title('Original Phase')
subplot(2,2,4)
plot(wavelength,phase)
title('Interpolated Phase')

M(:,1)=wavelength;
M(:,2)=intensity;
M(:,3)=phase;
dlmwrite('T:/Tobias/Daten/20140701_Chirped_Mirror_FROGS/001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A.bin.Speck_equidistant_data.csv',M,'delimiter',',','precision',8)

clear