pattern = '(?<file>\d+)_x=(?<x>[\d.]+)_y=(?<y>[\d.]+)_I=(?<I>[\d.]+)_W=(?<W>[\d.]+)';
folder = '\\gar-sv-home01/Tobias.Pleyer/Desktop/Master_log/Data/Daten/20140509_x-y-data-spectrum';

files = dir(folder);
L = length(files);

avgs = zeros(1,L);

colors = {'red','blue','green','black','magenta','cyan'};

figure(1)

chosen_files = [10 18 15 26 21 35 42 46];
legends = {'','','','','','','',''};


for i=1:8
    subplot(4,2,i)
   file = files(chosen_files(i)).name;
   if strcmp(file,'.') || strcmp(file,'..'), continue, end
   info = regexp(file, pattern, 'names');
   data = dlmread([folder '/' file],'\t',18, 0);
   wavl = data(:,1);
   spec = data(:,2);
   sampling_freq = 1e6/(wavl(2)-wavl(1));
   order         = 2;
   cut_off_freq  = 1e5;
   peak_to_peak_dB = 0.1;
   stopband_atten  = 20;
   normed_cutoff = 0.05;
   [B,A] = ellip(order,peak_to_peak_dB,stopband_atten,normed_cutoff,'low');
   filtered_spec = filtfilt(B,A,spec);
   weighted_avg = sum(wavl.*filtered_spec)/sum(filtered_spec);
   offset = (sum(filtered_spec(1:1000))/1000 + sum(filtered_spec(2200:end))/length(filtered_spec(2200:end)))/2;
   corrected_filtered_spec = filtered_spec - offset;
   M = max(corrected_filtered_spec);
   normalized_corrected_filtered_spec = corrected_filtered_spec / M;
   avgs(i) = weighted_avg;
   j = mod(i,6);
   plot(wavl,normalized_corrected_filtered_spec,'color',colors{j+1})
   s = sprintf('x=%s, y=%s',info.x,info.y);
   legend(s)
%    legends{i} = s;
end
% axis([950 1100 -0.1 1.1])
% xlabel('wavelength')
% ylabel('normed intensity')
% title('Overlay of all normalized spectra')
% legend(legends{1},legends{2},legends{3},legends{4},legends{5},legends{6},legends{7},legends{8})
% 
% figure(2)
% scatter([1:20 22:48],avgs(3:end),'+')
% axis([0 50 990 1040])
% xlabel('File number')
% ylabel('weighted average')
% title('scatter plot of the weighted spectra')