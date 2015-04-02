%##########################################################
%
% Show the FROG results for the chirped mirror compressor.
% Step3: Find the new spectra
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

E1 = dlmread('001_SHG_FROG_Compressor_Front_Mirror_Two_Bounces.bin.Ek.dat');
S1 = dlmread('001_SHG_FROG_Compressor_Front_Mirror_Two_Bounces.bin.Speck.dat');
S2 = dlmread('002_SHG_FROG_Compressor_Front_Mirror_Two_Bounces.bin.Speck.dat');
S1_1 = fftc(sqrt(fliplr(E1(:,2)')) .* exp(1i.*fliplr(E1(:,3)')));

figure()
hold on
plot(S1(:,2),'b-')
plot(S2(:,2),'g-')
plot(abs(S1_1).^2,'r-')
hold off