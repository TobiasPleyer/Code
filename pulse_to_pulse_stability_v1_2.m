%##########################################################################
%          This is a script to quantify the pulse to pulse
%          stability of a amplified laser.
%          
%          Version: 1.2
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
%##########################################################################

times = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_times_test.txt','w');
data = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_data_test.txt','w');

log_line = 'Time;Ampl';
content = fileread('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000.txt');
ends = regexp(content,log_line,'end');
fprintf(times,content(1:ends-9));
fprintf(data,content(ends-8:end));
fclose('all');

calc = false;
if calc
    % Start the calculations
    times = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_times_test.txt','r');
    data = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_data_test.txt','r');
    ptp = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_ptp2.txt','w');

    tmp = fgetl(times); % get rid of unnessecary header lines
    tmp = fgetl(times); %          ------""------
    tmp = fgetl(times); %          ------""------
    tmp = fgetl(data);  %          ------""------
    clear tmp

    time_curve = zeros(1,101);
    ampl_curve = zeros(1,101);
    for i=1:5000
        if mod(i,50)==0
            fprintf('%d\n',i)
        end
        log_line = fgetl(times);
        C = strsplit(log_line,';');
        log_value = str2double(C(end));
        for j=1:101
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
        p = polyfit(time,zoomed_ampl,10);
        time = zoomed_time(1)*1e9:0.005:zoomed_time(end)*1e9;
        dp = polyder(p);
        r = roots(dp); % find all roots of the derivative of the polynomial
        % We have to find the only reasonible x-value for us -> close to idx
        candidates = r(imag(r)==0);
        x = candidates(find_closest_idx(candidates,time_curve(idx)*1e9));
        theoretical_max = polyval(p,x);
        fprintf(ptp,'%f;%f\n',log_value,theoretical_max);
    end

    fclose('all');
end

result = dlmread('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_ptp.txt',';');
timestamps = result(:,1);
amplitudes = result(:,2);
differences = diff(amplitudes);
rms_ = rms(result(:,2));
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
