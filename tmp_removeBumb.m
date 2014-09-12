close all;

folder = '../Daten/FROGS 10kHz_2mJ/4mJ/';
fileBase = '001_AIR_FROG_20.6W_RTT=6.404us_Ip=12.0A.bin';

[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk] = compensation_loadData(folder,fileBase);

figure()
plotyy(l_Sk,I_Sk,l_Sk,p_Sk)

idx_low = find_closest_idx(1026.7,l_Sk);
idx_high = find_closest_idx(1028.5,l_Sk);

altered_p_Sk = p_Sk;
altered_p_Sk(idx_high:idx_low) = altered_p_Sk(idx_high:idx_low) + pi;

figure()
plotyy(l_Sk,I_Sk,l_Sk,altered_p_Sk)

altered_p_Sk2 = p_Sk;
L = abs(idx_high-idx_low) + 1;
factors = linspace(1,0.5,L);

altered_p_Sk2(idx_high:idx_low) = altered_p_Sk2(idx_high:idx_low) + factors'.*pi;

figure()
plotyy(l_Sk,I_Sk,l_Sk,altered_p_Sk2)

figure()
plot(l_Sk,altered_p_Sk2,'r')
hold on
plot(l_Sk,p_Sk)
hold off

dlmwrite('../Daten/FROGS 10kHz_2mJ/4mJ/complex_pulse_no_bumb_phase',[l_Sk,altered_p_Sk])
dlmwrite('../Daten/FROGS 10kHz_2mJ/4mJ/complex_pulse_no_bumb_intensity',[l_Sk,I_Sk])