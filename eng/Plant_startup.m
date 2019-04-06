[fd,~] = fopen('runs.dat','r');
if fd > 0
    tline = fgetl(fd);
    while ischar(tline)
        fprintf(1,'Processing: "%s"\n', tline);
        if exist(tline,'dir') == 7
            oldfolder = cd(tline);
            csv2mat;
            delete *.csv
            cd(oldfolder);
        end
        tline = fgetl(fd);
    end
    fclose(fd);
    delete runs.dat
end
clear fd msg oldfolder tline
ui_plant;
