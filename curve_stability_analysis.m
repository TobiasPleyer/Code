clear
set(0,'DefaultFigureWindowStyle','docked')

parent = 'Daten/';
filebase = '26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin';
filename = sprintf('%s%s.%s',parent,filebase,'Speck.dat');

% Preparing data to save the changed phase
special = '-wings-flattened';
figure_name = sprintf('%s%s','26-AIR-FROG-31W-RTT=7us-Ip=12A-Speck',special);
figure_file = sprintf('%s%s','Daten/Pictures/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_PICTURES/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck',special);

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

idxs = [1017.5 1018 1024 1027 1029.6 1031 1032.5 1034.5 1036 1039 1039.5];
[x_c,y_c,x_sep] = flatten_phase_curve(wavelength,phase,idxs);
bumbs= [1 1 0 0 0 0 1 1];
[x,y] = add_bumb(x_c,y_c,wavelength,phase,x_sep,bumbs,[],[],[],[]);

% [ax,h1,h2] = plotyy(wavelength,intensity,[wavelength,x],[phase,y]);
% legend([h1,h2],'Original','Modified','Spectrum')

% Show the difference between the Catmull-Rom ansatz and Matlab's smooth
% function
% figure(2)
% s_phase = smooth(phase);
% plot(wavelength,phase,'color','red','LineWidth',3);
% hold 'on'
% plot(x_c,y_c,'color','green','LineWidth',3,'LineStyle','--');
% plot(wavelength,s_phase,'color','blue','LineWidth',3,'LineStyle','-.');
% legend('Original','Catmull','Smooth')
% hold 'off'

new_phase    = zeros(size(wavelength));
smooth_phase = reshape(smooth(phase),1,256);


for i=1:length(wavelength)
   idx = find_closest_idx(x,wavelength(i));
   if idx==-1
      new_phase(i) = phase(i);
   else
      new_phase(i) = y(idx);
   end
end
% Only to use Matlab's smooth function
%new_phase = smooth_phase;
%new_phase = phase;
%r = (rand(size(new_phase)) - 0.5) .* 0.1;
%new_phase = new_phase + r;

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
hold 'off'

% For visual verification how the original and the modified phase differ
figure(3)
clf
plot(wavelength,phase-new_phase,'color','red')
axis([1015 1040 -2 2])
title(sprintf('Absolute difference between original and modified phase for\n%s',figure_file))
xlabel('wavelength in nm')
ylabel('absolute difference')
hold 'on'
L = length(idxs);
for i=1:L
    idx = find_closest_idx(wavelength,idxs(i));
    if idx==-1
            error('Index not in range!')
    end
    x_val = wavelength(idx);
    if (i>1 && i<L)
        line([x_val x_val],[-15 15],'LineStyle','--')
    end
end
hold 'off'

out(:,1) = wavelength;
out(:,2) = intensity;
out(:,3) = new_phase;       % the new phase
out(:,4) = phase-new_phase; % the effective compensated phase

figure(1)
saveas(gcf,sprintf('%s_%s',figure_file,'Fig1'),'jpg')
figure(3)
saveas(gcf,sprintf('%s_%s',figure_file,'Fig2'),'jpg')
parent = 'Daten/26_AIR_FROG_31W_RTT=7us_Ip=12A_Speck_CHANGED/';
output_file = sprintf('%s%s.%s%s.dat',parent,filebase,'Speck',special);
dlmwrite(output_file,out,'delimiter','\t','precision',6)








