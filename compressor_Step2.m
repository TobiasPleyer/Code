%##########################################################
%
% Show the FROG results for the chirped mirror compressor.
% Step2: Perform temporal flip to yield the same temporal
% pulses. Circular shifts to the spectrum
% Date: 1/8/2015
% Author: Tobias Pleyer
%
%##########################################################

N = 4; % number of consecutive measurements

file_counter = '%03d';
file_base = '_SHG_FROG_Compressor_Front_Mirror_Two_Bounces';
file_end = '.bin.Speck.dat';
folder_base = 'T:\LEX_measurements\hybrid data\';
folder_path = '20141218\Tobi\FROG\';

cd([folder_base folder_path]) % change to the working directory
colors = {'r-', 'b-', 'g-', 'c-', 'm-', 'k-'};
figure();

subplot(2,2,1)
hold on
leg = {}; % the legend object to be displayed
for n=1:N
    file_name = sprintf([file_counter file_base file_end],n);
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    if (n==1) || (n==4)
        plot(S(:,1),fftshift(S(:,2)),colors{n})
        leg{n} = sprintf([file_counter '+ fftshift'],n);
    else
        plot(S(:,1),S(:,2),colors{n})
        leg{n} = sprintf(file_counter,n);
    end
end
legend(leg,'Location','SouthEast')
title('Plot of the spectrum of the measured pulses','Fontsize',18)
xlabel('wavelength [nm]')
ylabel('a.u.')
hold off

subplot(2,2,2)
hold on
leg = {}; % the legend object to be displayed
for n=1:N
    file_name = sprintf([file_counter file_base file_end],n);
    try
        S = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    if (n==1) || (n==4)
        plot(S(:,1),(-1).*fftshift(S(:,3)),colors{n})
        leg{n} = sprintf([file_counter '+ fftshift + *(-1)'],n);
    else
        plot(S(:,1),S(:,3),colors{n})
        leg{n} = sprintf(file_counter,n);
    end
end
legend(leg,'Location','SouthEast')
title('Plot of the spectral phase of the measured pulses','Fontsize',18)
xlabel('wavelength [nm]')
ylabel('a.u.')
hold off

file_end = '.bin.Ek.dat';

subplot(2,2,3)
hold on
leg = {}; % the legend object to be displayed
for n=1:N
    file_name = sprintf([file_counter file_base file_end],n);
    try
        E = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    if (n==1) || (n==4)
        plot(E(:,1)',fliplr(E(:,2)'),colors{n})
        leg{n} = sprintf([file_counter '+ time flip'],n);
    else
        plot(E(:,1),E(:,2),colors{n})
        leg{n} = sprintf(file_counter,n);
    end
end
legend(leg,'Location','SouthEast')
title('Plot of the temporal form of the measured pulses','Fontsize',18)
xlabel('time [fs]')
ylabel('a.u.')
hold off

subplot(2,2,4)
hold on
leg = {}; % the legend object to be displayed
for n=1:N
    file_name = sprintf([file_counter file_base file_end],n);
    try
        E = dlmread(file_name);
    catch err
        fprintf('Error: File not found\n')
        continue
    end
    if (n==1) || (n==4)
        plot(E(:,1),(-1).*fftshift(E(:,3)),colors{n})
        leg{n} = sprintf([file_counter '+ fftshift + *(-1)'],n);
    else
        plot(E(:,1),E(:,3),colors{n})
        leg{n} = sprintf(file_counter,n);
    end
end
legend(leg,'Location','SouthEast')
title('Plot of the temporal phase of the measured pulses','Fontsize',18)
xlabel('time [fs]')
ylabel('a.u.')
hold off