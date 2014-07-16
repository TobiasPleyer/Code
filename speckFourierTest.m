Ek = dlmread('Daten\26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Ek.dat');
Sk = dlmread('Daten\26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck.dat');
% 
% wavelength = Sk(:,1);
% wavelength = reshape(wavelength,1,256);
% wavelength = fliplr(wavelength);
% 
% intensity = Sk(:,2);
% intensity = reshape(intensity,1,256);
% intensity = fliplr(intensity);
% 
% phase = Sk(:,3);
% phase = reshape(phase,1,256);
% phase = fliplr(phase);

Sk_cplx = sqrt(Sk(:,2)) .* exp(1i*Sk(:,3));

[t2,Ek2] = Speck_Fourier(Sk(:,1),Sk_cplx);

figure(1)
plot(Ek(:,1),Ek(:,2))
figure(2)
plot(t2,abs(Ek2),'color','red')

figure(3)
c = 299792458;
f = c./wavelength;
f0 = c/1030;
S = exp(-(f-f0).^2/(0.001*f0)^2);
S = S./max(S);
plot(wavelength,S)
figure(4)
[t,E] = Speck_Fourier(wavelength,S);
plot(t,abs(E))