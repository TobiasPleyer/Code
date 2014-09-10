function [t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]=compensation_loadData(fileBase, fileName)
    %% This is the part where all filenames and directories are defined

    if nargin<2
        error('Usage: compensation_loadData(fileBase, file1[, file2, file3, ...]\n')
    end
        filename_Et    = sprintf('%s%s.Ek.dat',fileBase,fileName);
        filename_Speck = sprintf('%s%s.Speck.dat',fileBase,fileName);
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