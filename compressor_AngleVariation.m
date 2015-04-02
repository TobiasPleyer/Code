filename = 'T:\Tobias\Chirped Mirror Compressor\Originals\Trubetskov_Mirror_Specification_Broad_Band.txt';

temp = dlmread(filename,'\t',1,0);
wavelength = temp(:,1);
GDD = temp(:,4);

before = 995:0.1:1009.9;
after = 1050.1:0.1:1065;
extended_wavelength = [before' ; wavelength ; after'];
extended_GDD = [zeros(150,1) ; GDD ; zeros(150,1)];

N = 8;
lower = -25;
upper = 25;
shift_array = round(linspace(lower,upper,N));
random_shift_array = datasample(shift_array,N);
Summe_linear = zeros(length(extended_GDD),1);
Summe_random = zeros(length(extended_GDD),1);

figure()
hold on
for i=1:N
    temp = circshift(extended_GDD,shift_array(i));
    plot(extended_wavelength,temp)
    Summe_linear = Summe_linear + temp;
end
hold off
title('The phase curves of the misaligned mirrors. The error was taken linear.','Fontsize',14)
xlabel('wavelength [nm]','Fontsize',14)
ylabel('phase [rad]','Fontsize',14')

figure()
hold on
for i=1:N
    temp = circshift(extended_GDD,random_shift_array(i));
    plot(extended_wavelength,temp)
    Summe_random = Summe_random + temp;
end
hold off
title('The phase curves of the misaligned mirrors. The error was taken random.','Fontsize',14)
xlabel('wavelength [nm]','Fontsize',14)
ylabel('phase [rad]','Fontsize',14)

figure()
plot(extended_wavelength,N.*extended_GDD)
hold on
plot(extended_wavelength,Summe_linear,'r-')
plot(extended_wavelength,Summe_random,'g-')
hold off
title('Difference between the unaltered phase curve in the case of linear and random alignment error','Fontsize',14)
legend('No Error','Linear Error','Random Error')
xlabel('wavelength [nm]','Fontsize',14)
ylabel('phase [rad]','Fontsize',14)

