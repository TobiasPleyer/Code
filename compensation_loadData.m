function [t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]=compensation_loadData()
    %% This is the part where all filenames and directories are defined

        parent         = '../Daten/Chirped Mirrors/Frogs/25.3W_RTT=6.404us_Ip=12.8A/';
        filebase       = '001_AIR_FROG_25.3W_RTT=6.404us_Ip=12.8A.bin';  
        filename_Et    = sprintf('%s%s.Ek.dat',parent,filebase);
        filename_Speck = sprintf('%s%s.Speck.dat',parent,filebase);

    %%


    %% Open the file and get the original data

        % Time based field   
        Et   = dlmread(filename_Et);
        t_Et = Et(:,1);
        I_Et = Et(:,2);
        p_Et = Et(:,3);

        % Wavelength based field
        Sk   = dlmread(filename_Speck);
        l_Sk = Sk(:,1);
        I_Sk = Sk(:,2);
        p_Sk = Sk(:,3);

    %%
end