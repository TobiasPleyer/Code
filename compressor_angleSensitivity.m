
str1 = 'T:\Tobias\Chirped Mirror Compressor\Analysis\002_SHG-FROG_Compressor_angle_check_M';
str2 = '.bin.Speck.dat';
w = {};
p = {};
Int = {};

for i=1:12
    name = [str1 int2str(i) str2];
    M = dlmread(name);
    w{i} = M(:,1);
    Int{i} = M(:,2);
    p{i} = M(:,3);
end

%The following lines are due to the issue that two files exhibited fourier
%shifts and/or sign problems. In order to keep the original data untouched
%I edited the files and saved them as .mat files. The unaltered data is
%substituted by this data.
load('T:\Tobias\Chirped Mirror Compressor\Analysis\M2')
load('T:\Tobias\Chirped Mirror Compressor\Analysis\M7')
w{2} = M2(:,1);
Int{2} = M2(:,2);
p{2} = M2(:,3);
w{7} = M7(:,1);
Int{7} = M7(:,2);
p{7} = M7(:,3);

common_wavelength = 1000:0.1:1070;

for i=1:12
    p{i} = interp1(w{i},p{i},common_wavelength);
    Int{i} = interp1(w{i},Int{i},common_wavelength);
end

for i=2:12
    p{i} = compressor_alignCurve(common_wavelength,p{1},p{i});
end

colors = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};
plot(common_wavelength,p{i},'Color',colors{mod(1,6)+1},'LineWidth',4)
hold on
for i=2:12
    plot(common_wavelength,p{i},'Color',colors{mod(i,6)+1})
end
%title(sprintf('Effect of different angles of incidence. The wheel of the\nmirror mount was varied between 0 and 1000 degrees.'),'FontSize',18)
xlabel('wavelength [nm]','FontSize',18)
ylabel('phase [rad]','FontSize',18)
xlim([1010 1050])
set(gca,'FontSize',16)
hold off