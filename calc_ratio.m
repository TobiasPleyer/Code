%##########################################################################
%          This is a script to quantify the power output from our
%          laser diodes.
%          
%          Version: 1.0
%
%          Output: A semicolon (';') seperated txt file for each
%                  calibration file provided.
%          Format: Timestamp;Amplitude
%          
%          Values of zero or below are linearly fitted using Matlab's built
%          in interp function.
%
%          HISTORY
%          ----------------------------------------------------------------
%          v1.0: First prototype of the program (2014/06/16)
%
%##########################################################################

data1 = dlmread('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_1_no_head.txt','\t',1,0);
data2 = dlmread('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_2_no_head.txt','\t',1,0);

timestamp1 = data1(:,1);
C_A1 = data1(:,2);
C_B1 = data1(:,3);
timestamp2 = data2(:,1);
C_A2 = data2(:,2);
C_B2 = data2(:,3);

clear data1; clear data2
% fill in missing data
l = length(C_A1);
first_zero = 1;
next_nonzero = 1;
while first_zero<=l
    if C_A1(first_zero)==0
        next_nonzero = first_zero;
        while C_A1(next_nonzero)==0 && next_nonzero<l
            next_nonzero = next_nonzero + 1;
        end
        mid = (C_A1(next_nonzero)+C_A1(first_zero-1))/2;
        C_A1(first_zero:next_nonzero-1) = mid;
        first_zero = next_nonzero + 1;
    else
        first_zero = first_zero + 1;
    end
end
l = length(C_B1);
first_zero = 1;
next_nonzero = 1;
while first_zero<=l
    if C_B1(first_zero)==0
        next_nonzero = first_zero;
        while C_B1(next_nonzero)==0 && next_nonzero<l
            next_nonzero = next_nonzero + 1;
        end
        mid = (C_B1(next_nonzero)+C_B1(first_zero-1))/2;
        C_B1(first_zero:next_nonzero-1) = mid;
        first_zero = next_nonzero + 1;
    else
        first_zero = first_zero + 1;
    end
end
l = length(C_A2);
first_zero = 1;
next_nonzero = 1;
while first_zero<=l
    if C_A2(first_zero)==0
        next_nonzero = first_zero;
        while C_A2(next_nonzero)==0 && next_nonzero<l
            next_nonzero = next_nonzero + 1;
        end
        mid = (C_A2(next_nonzero)+C_A2(first_zero-1))/2;
        C_A2(first_zero:next_nonzero-1) = mid;
        first_zero = next_nonzero + 1;
    else
        first_zero = first_zero + 1;
    end
end
l = length(C_B2);
first_zero = 1;
next_nonzero = 1;
while first_zero<=l
    if C_B2(first_zero)==0
        next_nonzero = first_zero;
        while C_B2(next_nonzero)==0 && next_nonzero<l
            next_nonzero = next_nonzero + 1;
        end
        mid = (C_B2(next_nonzero)+C_B2(first_zero-1))/2;
        C_B2(first_zero:next_nonzero-1) = mid;
        first_zero = next_nonzero + 1;
    else
        first_zero = first_zero + 1;
    end
end

data1(:,1) = timestamp1;
data1(:,2) = C_A1;
data1(:,3) = C_B1;
data1(:,4) = C_A1./C_B1;
data2(:,1) = timestamp2;
data2(:,2) = C_A2;
data2(:,3) = C_B2;
data2(:,4) = C_A2./C_B2;
dlmwrite('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_1_ratio.txt',data1,'delimiter','\t');
dlmwrite('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_2_ratio.txt',data2,'delimiter','\t');
clear data1; clear data2