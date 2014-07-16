clear
set(0,'DefaultFigureWindowStyle','docked')

parent = 'Daten/';
filebase = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
filename = sprintf('%s%s.%s',parent,filebase,'Speck.dat');

% Preparing data to save the changed phase
special = 'adjust-height-10268-10275';
figure_name = sprintf('%s%s','26-AIR-FROG-31W-RTT=7us-Ip=12A-Speck',special);
figure_file = sprintf('%s-%s','Daten/Pictures/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_PICTURES/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck',special);

fhandle = dlmread(filename);

wavelength = fhandle(:,1);
wavelength = reshape(wavelength,1,256);
wavelength = fliplr(wavelength);

intensity = fhandle(:,2);
intensity = reshape(intensity,1,256);
intensity = fliplr(intensity);

phase = fhandle(:,3);
phase = reshape(phase,1,256);
phase = fliplr(phase);

% The brains of this whole script
w1 = 1026.8;
w2 = 1027.5;
new_phase = adjust_height(wavelength,phase,w1,w2);
new_phase = reshape(new_phase,1,256);

figure(1)
plot(wavelength,phase,'color','red')
title(sprintf('View of the whole phase curve and the smoothened area for\n%s',figure_name))
xlabel('wavelength in nm')
ylabel('phase value')
axis([1018 1040 -15 15])
hold 'on'
plot(wavelength,new_phase)
temp = (intensity.*28)-15;
plot(wavelength,temp,'color','green')
legend('original','modified','intensity','Location','NorthWest')
line([w1 w1],[-15 15],'LineStyle','--')
line([w2 w2],[-15 15],'LineStyle','--')
hold 'off'

% For visual verification how the original and the modified phase differ
figure(3)
clf
plot(wavelength,phase-new_phase,'color','blue')
axis([1015 1040 -5 5])
title(sprintf('Absolute difference between original and modified phase for\n%s',figure_file))
xlabel('wavelength in nm')
ylabel('absolute difference')


out(:,1) = wavelength;
out(:,2) = intensity;
out(:,3) = new_phase;       % the new phase
out(:,4) = phase-new_phase; % the effective compensated phase

figure(1)
saveas(gcf,sprintf('%s-%s',figure_file,'Fig1'),'jpg')
figure(3)
saveas(gcf,sprintf('%s-%s',figure_file,'Fig2'),'jpg')
parent = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/';
output_file = sprintf('%s%s.%s-%s.dat',parent,filebase,'Speck',special);
dlmwrite(output_file,out,'delimiter','\t','precision',6)