% This script calculates the portion of the pulse's energy which is not
% located in the central peak of the pulse.
% This portion is desired to be as small as possible, where the respective
% portion of the Fourier limit imposes the physically achievable limit.

clear

dir_main = '\\gar-sv-home01\Tobias.Pleyer\Desktop\Matlab Moritz\FROG\';
dir_spec = 'Spectrum/';
dir_time = 'Time/';
dir_data = 'Data/';
filename_base = '09_AIR_FROG_18.0W_RTT=7.415us_Ip=10.0A.bin';

phase_factor = 1;
tF1 = -220;
tF2 = -tF1;
t1 = -270;                                                                 % Manual data for the boundaries of the main peak.
t2 = 200;                                                                  % Everything outside these borders is part of the pulse wings
c = 299792458; %m/s

filename_Ek = sprintf('%s%s.Ek.dat',dir_main, filename_base);
Ek = dlmread(filename_Ek);

t = Ek(:,1);
I = Ek(:,2);

ind1 = find(t >= t1, 1, 'first');
ind2 = find(t >= t2, 1, 'first');

I0 = I/trapz(t, I);                                                       % Calculate the estimated trapezoidal integral of the I with respect to t
I_side = trapz(t(ind1:ind2), I0(ind1:ind2));

figure(1)
plot(t, I0)
hold on
plot([t(ind1), t(ind1)],[0, max(I0)],'r')
plot([t(ind2), t(ind2)],[0, max(I0)],'r')
hold off
s = sprintf('Energy of Peak %.0f%%',I_side*100);
title(s);
ylabel('Intensity (a.u.)');
xlabel('Time (fs)');


filename_Speck = sprintf('%s%s.Speck.dat',dir_main, filename_base);
Speck = dlmread(filename_Speck);
f = c./(Speck(:,1).*1e-9);
omega = 2.*pi.*f;
T = 1/(f(2)-f(1));
t = T/2.*linspace(-1,1,size(f,1)).*1e15;
t=t';

Sw = (Speck(:,1).*1e-9).^2./(2*pi*c).*Speck(:,2);
Sw = Sw./max(Sw);
Phiw = Speck(:,3);
Phiw = Phiw.*phase_factor;
A_1f = sqrt(Sw).*exp(-1i.*Phiw);
%Get Fourier limit
EkFourier=fftshift(fft(ifftshift(sqrt(Sw))));
IFourier = abs(EkFourier/max(abs(EkFourier))).^2;

indF1 = find(t >= tF1, 1, 'first');
indF2 = find(t >= tF2, 1, 'first');

IFourier = IFourier/trapz(t, IFourier);
I_sideF = trapz(t(indF1:indF2), IFourier(indF1:indF2));

figure(2)

plot(t, IFourier)
hold on
plot([t(indF1), t(indF1)],[0, max(IFourier)],'r')
plot([t(indF2), t(indF2)],[0, max(IFourier)],'r')
hold off
s = sprintf('Energy of Peak (Fourier limit) %.0f%%',I_sideF*100);
title(s);
ylabel('Intensity (a.u.)');
xlabel('Time (fs)');
