function [D1,D2,D3]=compensation_calcDevs(p_Sk,w_Sk,B,A)
    D1 = diff(p_Sk) ./ (w_Sk(1)-w_Sk(2));
    % Immediately apply our lowpass filter
    %D1 = filtfilt(B,A,D1);
    D2 = diff(D1) ./ (w_Sk(1)-w_Sk(2));
    D2 = filtfilt(B,A,D2);
    D3 = diff(D2) ./ (w_Sk(1)-w_Sk(2));
end