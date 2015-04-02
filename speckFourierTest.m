Ek = dlmread('001_SHG_FROG_Compressor_Front_Mirror_Two_Bounces.bin.Ek.dat');
Sk = dlmread('001_SHG_FROG_Compressor_Front_Mirror_Two_Bounces.bin.Speck.dat');

wavelength = Sk(:,1);
wavelength = reshape(wavelength,1,256);
wavelength = fliplr(wavelength);

intensity = Sk(:,2);
intensity = reshape(intensity,1,256);
intensity = fliplr(intensity);

phase = Sk(:,3);
phase = reshape(phase,1,256);
phase = fliplr(phase);

Sk_cplx = sqrt(Sk(:,2)) .* exp(1i*Sk(:,3));

[t2,Ek2] = Speck_Fourier(Sk(:,1),Sk_cplx);
Sk2 = fftc(sqrt(Ek(:,2)).*exp(1i.*Ek(:,3))); %fftc is a function from Rick Trebino's code
                                       %and has to be visible in order for
                                       %this to work

% figure()
% plot(Ek(:,1),Ek(:,2))
% figure()
% plot(t2,abs(Ek2),'color','red')
figure()
plot(wavelength,intensity,'b-')
hold on
plot(wavelength,abs(Sk2).^2,'r-')
hold off

% figure()
% c = 299792458;
% f = c./wavelength;
% f0 = c/1030;
% S = exp(-(f-f0).^2/(0.001*f0)^2);
% S = S./max(S);
% plot(wavelength,S)
% figure()
% [t,E] = Speck_Fourier(wavelength,S);
% plot(t,abs(E))