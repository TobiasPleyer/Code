%##########################################################################
%          This is a script to quantify the pulse to pulse
%          stability of a amplified laser.
%          
%          Version: 2.0
%
%          Output: A semicolon (';') seperated txt file with two columns
%                  Timestamp;Amplitude
%          
%          The amplitude values are obtained by fitting a polynomial to the
%          curve.
%
%          HISTORY
%          ----------------------------------------------------------------
%          v1.0: First prototype of the program. (2014/06/18)
%          v1.1: Removed the need to split the original file. (2014/06/18)
%          v1.2: Removed all commentes and unnessecary code. (2014/06/18)
%          v1.3: Removed the filenames in favor of a filename variable.
%                Added more statistics (std, mean, etc.) (2014/06/23)
%          v2.0: Made the program smarter and added explanations. I also
%                reduced the degree of the polynomial. (2014/06/26)
%##########################################################################

clear; warning off;

plt = false;
calc = true;

file_base_name = '../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000';

times = fopen([file_base_name '_times.txt'],'w');
data = fopen([file_base_name '_data.txt'],'w');

log_line = 'Time;Ampl';
content = fileread([file_base_name '.txt']);
ends = regexp(content,log_line,'end');
cut = 0;
% ASCII INFO
% 10 = new line
% 13 = carriage return
% 32 = white space
while true
    if uint8(content(ends-9-cut))==10 || uint8(content(ends-9-cut))==13 || uint8(content(ends-9-cut))==32
        cut = cut+1;
    else
        break
    end
end
fprintf(times,content(1:ends-9-cut));
cut = 0;
while true
    if uint8(content(end-cut))==10 || uint8(content(end-cut))==13 || uint8(content(end-cut))==32
        cut = cut+1;
    else
        break
    end
end
fprintf(data,content(ends-8:end-cut));
fclose('all');

if calc
    % Start the calculations
    times = fopen([file_base_name '_times.txt'],'r');
    data = fopen([file_base_name '_data.txt'],'r');

    tmp = fgetl(times); % get rid of unnessecary header lines
    tmp = fgetl(times); % read header information
    idxs = strfind(tmp,';');
    segments = str2num(tmp(idxs(1)+1:idxs(2)-1));
    size = str2num(tmp(idxs(3)+1:end));
    tmp = fgetl(times); % get rid of unnessecary header lines
    tmp = fgetl(data);  %          ------""------
    clear tmp

    time_curve  = zeros(1,size);
    ampl_curve  = zeros(1,size);
    theo_values = zeros(1,segments);
    log_values  = zeros(1,segments);
    for i=1:segments
        if mod(i,250)==0
            fprintf('%d\n',i)
        end
        log_line = fgetl(times);
        C = strsplit(log_line,';');
        log_values(i) = str2double(C(end));
        for j=1:size
            value_line = fgetl(data);
            C = strsplit(value_line,';');
            time_curve(j) = str2double(C(1));
            ampl_curve(j) = str2double(C(2));
        end
        % Here starts the main logic. Try to fit a parabola to measured data.
        [m,idx] = max(ampl_curve);
        zoomed_ampl = ampl_curve(ampl_curve./m>0.2); % get rid of the unimportant data in order to properly fit the curve
        zoomed_time = time_curve(ampl_curve./m>0.2);
        time = zoomed_time*1e9;
        p = polyfit(time,zoomed_ampl,6);
        time = zoomed_time(1)*1e9:0.005:zoomed_time(end)*1e9;
        dp = polyder(p);
        r = roots(dp); % find all roots of the derivative of the polynomial
        % We have to find the only reasonible x-value for us -> close to idx
        candidates = r(imag(r)==0);
        x = candidates(find_closest_idx(candidates,time_curve(idx)*1e9));
        theo_values(i) = polyval(p,x);
    end
    ptp = fopen([file_base_name '_ptp2.txt'],'w');
    mean_ = mean(theo_values);
    std_ = std(theo_values);
    fprintf(ptp,'# mean: %f\n',mean_);
    fprintf(ptp,'# std: %f\n',std_);
    fprintf(ptp,'# Segments: %d\n',segments);
    fprintf(ptp,'# %s;%s\n','Timestamp','Peak');
    for i=1:segments
        fprintf(ptp,'%f;%f\n',log_values(i),theo_values(i));
    end

end
fclose('all');

if plt
    result = dlmread([file_base_name '_ptp.txt'],';');
    timestamps = result(:,1);
    amplitudes = result(:,2);
    differences = diff(amplitudes);
    percentages = differences./amplitudes(1:end-1)*100;
    figure(1)
    plot(timestamps,amplitudes)
    text(0.5,0.71,'rms: 0.6720')
    figure(2)
    plot(timestamps(1:end-1),differences)
    title('Difference from one pulse to its neighbour')
    ylabel('difference')
    xlabel('Timestamp')
    figure(3)
    plot(timestamps(1:end-1),percentages)
    title('Difference from one pulse to its neighbour in percent')
    ylabel('difference in percent compared to first pulse')
    xlabel('Timestamp')
end

warning on