%##########################################################################
%          This is a script to quantify the power output from our
%          laser diodes.
%          
%          Version: 2.0
%
%          Output: A semicolon (';') seperated txt file for each
%                  calibration file provided.
%          Format: Timestamp;Amplitude
%          
%          Values of zero or below are linearly fitted using Matlab's built
%          in interp1 function.
%
%          HISTORY
%          ----------------------------------------------------------------
%          v1.0: First prototype of the program (2014/06/16)
%          
%          v2.0: Introduced a major simplification of the implementation.
%                (2014/06/23)
%
%##########################################################################

data1 = dlmread('t:\LEX_measurements\hybrid data\20140712\CalibrationRegen_nohead.txt','\t',1,0);
%data2 = dlmread('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_2_no_head.txt','\t',1,0);

timestamp1 = data1(:,1);
C_A1 = data1(:,2);
C_B1 = data1(:,3);
%timestamp2 = data2(:,1);
%C_A2 = data2(:,2);
%C_B2 = data2(:,3);

data1(:,2) = interp1(timestamp1(C_A1 > 0),C_A1(C_A1 > 0),timestamp1,'linear','extrap');
data1(:,3) = interp1(timestamp1(C_B1 > 0),C_B1(C_B1 > 0),timestamp1,'linear','extrap');
%data2(:,2) = interp1(timestamp2(C_A2 > 0),C_A2(C_A2 > 0),timestamp2,'linear','extrap');
%data2(:,3) = interp1(timestamp2(C_B2 > 0),C_B2(C_B2 > 0),timestamp2,'linear','extrap');

data1(:,4) = data1(:,2)./data1(:,3);
%data2(:,4) = data2(:,2)./data2(:,3);
dlmwrite('t:\LEX_measurements\hybrid data\20140712\CalibrationRegen_nohead_ratio.txt',data1,'delimiter',';');
%dlmwrite('\\gar-sv-home01\Tobias.Pleyer\Desktop\Master_log\Data\Daten\20140616_Diodenmessungen\Calibration_2_ratio2.txt',data2,'delimiter','\t');
clear data1; 
%clear data2