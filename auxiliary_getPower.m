function [power]=auxiliary_getPower(filename)
    fid = fopen(filename);
    for i=1:18
        line = fgetl(fid);
    end
    fclose(fid);
    power = str2double(line(9:15));
end