%##########################################################################
%          This is a script to quantify the pulse to pulse
%          stability of a amplified laser.
%          
%          Version: 1.1
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
%##########################################################################

%file = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000.txt','r');
times = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_times_test.txt','w');
data = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_data_test.txt','w');
% first_header_lines  = 0;
% second_header_lines = 0;
% in_first_part       = false;
% in_second_part      = false;
% first_header_done   = false;
% second_header_done  = false;
% headers_done        = false;
% data_start          = false;
% log_lines           = 0;
% points_per_log      = 0;
% regular_expression  = '^#\d+';
% % Start to seperate the oscilloscope file into a log file and a data file
% while ~feof(file)
%     line = fgets(file);
%     if isempty(regexp(line,regular_expression,'once'))
%         if headers_done
%             error('Something is wrong in the file.')
%         elseif ~first_header_done
%             first_header_lines = first_header_lines + 1;
%             fprintf(times,line);
%         elseif ~second_header_done
%             in_first_part = false;
%             second_header_lines = second_header_lines + 1;
%             fprintf(data,line);
%         else
%             error('Something went wrong!')
%         end
%     elseif in_first_part
%         fprintf(times,line);
%     elseif in_second_part
%         fprintf(data,line);
%     elseif ~first_header_done
%         first_header_done = true;
%         in_first_part = true;
%         fprintf(times,line);
%     elseif ~second_header_done
%         second_header_done = true;
%         in_second_part = true;
%         headers_done = true;
%         fprintf(data,line);
%     else
%         error('Something is wrong in the file.')
%     end
% end
% MUCH FASTER!!!
log_line = 'Time;Ampl';
content = fileread('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000.txt');
ends = regexp(content,log_line,'end');
fprintf(times,content(1:ends-9));
fprintf(data,content(ends-8:end));
fclose('all');

% Start the calculations
% times = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_times.txt','w');
% data = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_data.txt','w');
% %ptp = fopen('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_ptp2.txt','w');
% 
% tmp = fgetl(times); % get rid of unnessecary header lines
% tmp = fgetl(times); %          ------""------
% tmp = fgetl(times); %          ------""------
% tmp = fgetl(data);  %          ------""------
% clear tmp
% 
% time_curve = zeros(1,101);
% ampl_curve = zeros(1,101);
% figure(1)
% for i=1:log_lines
%     if mod(i,50)==0
%         fprintf('%d\n',i)
%     end
%     log_line = fgetl(times);
%     C = strsplit(log_line,';');
%     log_value = str2double(C(end));
%     for j=1:points_per_log
%         value_line = fgetl(data);
%         C = strsplit(value_line,';');
%         time_curve(j) = str2double(C(1));
%         ampl_curve(j) = str2double(C(2));
%     end
%     % Here starts the main logic. Try to fit a parabola to measured data.
%     [m,idx] = max(ampl_curve);
%     zoomed_ampl = ampl_curve(ampl_curve./m>0.2); % get rid of the unimportant data in order to properly fit the curve
%     zoomed_time = time_curve(ampl_curve./m>0.2);
%     time = zoomed_time*1e9;
%     p = polyfit(time,zoomed_ampl,10);
%     time = zoomed_time(1)*1e9:0.005:zoomed_time(end)*1e9;
%     dp = polyder(p);
%     r = roots(dp); % find all roots of the derivative of the polynomial
%     % We have to find the only reasonible x-value for us -> close to idx
%     candidates = r(imag(r)==0);
%     x = candidates(find_closest_idx(candidates,time_curve(idx)*1e9));
%     theoretical_max = polyval(p,x);
%     curve = p(1)*time.^10 + ...
%             p(2)*time.^9  + ...
%             p(3)*time.^8  + ...
%             p(4)*time.^7  + ...
%             p(5)*time.^6  + ...
%             p(6)*time.^5  + ...
%             p(7)*time.^4  + ...
%             p(8)*time.^3  + ...
%             p(9)*time.^2  + ...
%             p(10)*time    + ...
%             p(11);
%     %fprintf(ptp,'%f;%f\n',log_value,theoretical_max);
% %     if mod(i,50)==0
% %         subplot(3,2,i/50)
% %         hold on;plot(zoomed_time,zoomed_ampl,'r');plot(time*1e-9,curve,'g');hold off;
% %     end
% end
% 
% fclose('all');
% 
% % result = dlmread('../Daten/20140617_Puls_zu_Puls/C3Stability_at_25.5W_5.1mJ_5kHz_200000_ptp.txt',';');
% % plot(result(:,1),result(:,2))