%% Preamble

% Use this script for general workup of chamber blanks

%% SETTINGS AND FOLDER INITIALIZATION
% Settings for a run are located in a user-defined config.ini file. Please
% see an example config.ini file for more details.

% Load configuration file containing settings
s = ini2struct('config.ini');
fields = fieldnames(s); % List field names in the settings structure
for i = 1:numel(fields) 
    s.(fields{i}) = str2double(s.(fields{i})); % Convert character strings to numerical values
end

s.plant_run_date = num2str(s.plant_run_date);
s.hcho_run_date  = num2str(s.hcho_run_date);
Plant_dir = ['D:\Plant\RAW\',s.plant_run_date,'\'];
HCHO_dir  = ['D:\Data\HCHO\RAW\',s.hcho_run_date,'\'];
addpath(Plant_dir)
addpath(HCHO_dir)

% Load data files
p = load('planteng_1.mat'); %Loads all variables from plant computer to the 'p' structure
load(['D:\Data\HCHO\RAW\',s.hcho_run_date,'\FILIF_ProcessedHCHO.mat']);  %% FIX THIS TO AUTOMATICALLY LOAD

%% Time

% Convert Plant Posixtime to Datetime
p.plant_datetime = datetime(p.Tplanteng_1,'ConvertFrom','posixtime');

% Convert to local time if specified in config.ini file
if s.local_time_convert
    p.plant_datetime = p.plant_datetime - hours(s.time_adjust);
    FILIF.datetime = FILIF.datetime - hours(s.time_adjust);
end

% Create posixtime array for FILIF and correct plant posixtime array for
% possible change due to local time conversion
FILIF.posixtime = posixtime(FILIF.datetime);
p.Tplanteng_1 = posixtime(p.plant_datetime);

%% TIME AVERAGE HCHO DATA AND INTERPOLATE ONTO PLANT TIME BASIS

[FILIF.time_1Hz, FILIF.hcho_1Hz] = binavg_FILIF(FILIF.posixtime,FILIF.hcho,1);
p.hcho = interp1(FILIF.time_1Hz,FILIF.hcho_1Hz,p.Tplanteng_1);


%% REMOVE DATA BY GRAPHICAL SELECTION (OPTIONAL)
% Manually remove points by using the brush tool provided in the Matlab
% figure options. Allows for easy removal of data when no readily-definable
% criterion exists to remove the data. Do not use this ability to
% 'cherry-pick' the data!

file_exist_check = exist(fullfile(Plant_dir,'RemovedPoints.mat'), 'file');

if file_exist_check == 2
    load('RemovedPoints.mat')
    
    for i = 1:length(p.Tplanteng_1)
        if ptsRemoved(i)
            p.CO2_ppm(i)    = NaN;
            p.H2O_ppth(i)   = NaN;
            p.hcho(i)       = NaN;
            p.SHT31_RH(i)   = NaN;
            p.SHT31_Temp(i) = NaN;
            p.TS_Temp(i)    = NaN;
        end
    end 
    loaded_ptsRemoved = ptsRemoved;
end

if s.graphical_removal
    
    figure
    h = pan;
    h.Motion = 'horizontal';
    h.Enable = 'on';
    ax1 = subplot(3,2,1);
    hLines = plot(p.plant_datetime,p.CO2_ppm);

    ax2 = subplot(3,2,2);
    plot(p.plant_datetime,p.H2O_ppth)
    
    ax3 = subplot(3,2,3);
    plot(p.plant_datetime,p.hcho)
    
    ax4 = subplot(3,2,4);
    plot(p.plant_datetime,p.SHT31_RH)
    
    ax5 = subplot(3,2,5);
    plot(p.plant_datetime,p.SHT31_Temp)
    
    ax6 = subplot(3,2,6);
    plot(p.plant_datetime,p.TS_Temp)

    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')

    % Start brushing mode and wait for user to hit "Enter" when done
    brush on
    disp('Hit Enter in command window when done brushing')
    pause

    % Loop through each graphics object
    for k = 1:numel(hLines)
        % Check that the property is valid for that type of object
        % Also check if any points in that object are selected
        if isprop(hLines(k),'BrushData') && any(hLines(k).BrushData)
            % Output the selected data to the base workspace with assigned name
            ptsRemoved = logical(hLines(k).BrushData.');
        end
    end  
end

% This step makes sure that the loaded ptsRemoved isn't overwritten by the
% new points being removed
if file_exist_check == 2
    ptsRemoved = ptsRemoved|loaded_ptsRemoved; % logical OR
end

if exist('ptsRemoved','var') == 1
    
    save(fullfile(Plant_dir,'RemovedPoints.mat'),'ptsRemoved');

    for i = 1:length(p.Tplanteng_1)
        if ptsRemoved(i)
            p.CO2_ppm(i)    = NaN;
            p.H2O_ppth(i)   = NaN;
            p.hcho(i)       = NaN;
            p.SHT31_RH(i)   = NaN;
            p.SHT31_Temp(i) = NaN;
            p.TS_Temp(i)    = NaN;
        end
    end
end

clear('ax1','ax2','ax3','ax4','ax5','ax6','h','hLines','file_exist_check')


%% Parsing of Bypass vs Chamber
% Flag 51 - Sampling Bypass
% Flag 50 - Sampling Chamber

BypassIndices = find(p.Flag==51);
ChamberIndices = find(p.Flag==50);

[BypassIndices, BypassIndices_rm]  = RemoveFirstPoints(BypassIndices, 0);
[ChamberIndices, ChamberIndices_rm]  = RemoveFirstPoints(ChamberIndices, 0);

% Separate bypass data from chamber data. All times in p are at 1 Hz.

Bypass_datetime_1Hz  = p.plant_datetime(BypassIndices);
Bypass_posixtime_1Hz = p.Tplanteng_1(BypassIndices);
Bypass_CO2_1Hz       = p.CO2_ppm(BypassIndices);
Bypass_H2O_1Hz       = p.H2O_ppth(BypassIndices);
Bypass_HCHO_1Hz      = p.hcho(BypassIndices);

Chamber_datetime_1Hz  = p.plant_datetime(ChamberIndices);
Chamber_posixtime_1Hz = p.Tplanteng_1(ChamberIndices);
Chamber_CO2_1Hz       = p.CO2_ppm(ChamberIndices);
Chamber_H2O_1Hz       = p.H2O_ppth(ChamberIndices);
Chamber_HCHO_1Hz      = p.hcho(ChamberIndices);

% Averaging
time_avg = 10; % Averaging time for data in seconds

if s.engplot
    figure,plot(Bypass_datetime_1Hz,Bypass_HCHO_1Hz,'.','MarkerSize',20)
    hold on
    plot(Chamber_datetime_1Hz,Chamber_HCHO_1Hz,'.','MarkerSize',20)
end

%% CHAMBER HCHO

[Chamber_posixtime, Chamber_HCHO] = binavg_plant(Chamber_posixtime_1Hz, Chamber_HCHO_1Hz,time_avg);

% Remove NaNs
temp_matrix = [Chamber_posixtime' Chamber_HCHO];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Chamber_posixtime = temp_matrix(:,1);
Chamber_HCHO = temp_matrix(:,2);
Chamber_datetime = datetime(Chamber_posixtime,'ConvertFrom','posixtime');

% Chunk the Newly-Averaged Chamber Data to Prepare for Outlier Removal
a = size(Chamber_datetime,1);
transitions = find(diff(Chamber_datetime) > minutes(4));
transitions = [transitions; a];
chamber_chunks = ones(size(transitions,1),2);
for i = 1:size(transitions,1)
    chamber_chunks(i,2) = transitions(i);
end
for i = 1:size(transitions,1)-1
    chamber_chunks(i+1,1) = transitions(i)+1;
end

% OUTLIER REMOVAL for Chamber Measurements
% Remove outliers more than three local scaled MAD from the local median
% Note: movmean (removing outliers more than three local standard deviations
% from the local mean) is too strongly affected by outliers

figure,plot(Chamber_datetime,Chamber_HCHO,'.','MarkerSize',15)
hold on
outlier_logical_array = [];
for i = 1:size(chamber_chunks,1)
    temp = isoutlier(Chamber_HCHO(chamber_chunks(i,1):chamber_chunks(i,2)),'gesd');
    outlier_logical_array = [outlier_logical_array; temp];
end
for j=1:length(outlier_logical_array)
    if outlier_logical_array(j)==1
        Chamber_HCHO(j) = NaN;
    end
end
plot(Chamber_datetime,Chamber_HCHO,'.','MarkerSize',15)
title('Chamber Outlier Removal')

% Remove NaNs
Chamber_posixtime = posixtime(Chamber_datetime);
temp_matrix = [Chamber_posixtime Chamber_HCHO];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Chamber_posixtime = temp_matrix(:,1);
Chamber_HCHO = temp_matrix(:,2);
Chamber_datetime = datetime(Chamber_posixtime,'ConvertFrom','posixtime');

% Interpolate Flag2 (representing equilibrated steps) onto Chamber_datetime
% 30 represents times when steps are equilibrated; defined prior to experiment
Chamber_Flag2 = interp1(p.plant_datetime,p.Flag2,Chamber_datetime);
Chamber_Equilibrated_Indices = find(Chamber_Flag2==30);

% Chunk equilibrated indices to break into their respective steps
Chamber_Equilibrated_Chunks = chunker(Chamber_Equilibrated_Indices);

% Show Flag 2 on Data
if s.engplot
    figure,plot(Chamber_datetime,Chamber_HCHO)
    hold on
    plot(Chamber_datetime,Chamber_Flag2,'.','MarkerSize',20)
    title('Chamber HCHO with Flag2')
end

Chamber_HCHO_avg = [];
Chamber_HCHO_avg_std = [];

for i = 1:size(Chamber_Equilibrated_Chunks,1)       
    Chamber_HCHO_avg(i)      = mean(Chamber_HCHO(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2)));
    Chamber_HCHO_avg_std(i)  = std(Chamber_HCHO(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2)))/sqrt(size(Chamber_HCHO(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2)),1));
end

% Create groups based on equilibrated indices
chamber_groups = [];
for i = 1:size(Chamber_Equilibrated_Chunks,1)
    temp = size(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2),2);
    v = zeros(temp, 1);
    v(:) = i;
    chamber_groups = [chamber_groups; v];
end

% Pull out equilibrated chamber values
chamber_HCHO_pts = [];
for i = 1:size(Chamber_Equilibrated_Chunks,1)
    temp = Chamber_HCHO(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2));
    chamber_HCHO_pts = [chamber_HCHO_pts; temp];
end

%% BYPASS HCHO
% Obtain averaged bypass HCHO for each step of an experimental run with
% corresponding standard deviation of the mean

[Bypass_posixtime, Bypass_HCHO] = binavg_plant(Bypass_posixtime_1Hz, Bypass_HCHO_1Hz,time_avg);

% Remove NaNs
temp_matrix = [Bypass_posixtime' Bypass_HCHO];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Bypass_posixtime = temp_matrix(:,1);
Bypass_HCHO = temp_matrix(:,2);
Bypass_datetime = datetime(Bypass_posixtime,'ConvertFrom','posixtime');

% Chunk the Newly-Averaged Bypass Data to Prepare for Outlier Removal
a = size(Bypass_datetime,1);
transitions = find(diff(Bypass_datetime) > minutes(10));
transitions = [transitions; a];
bypass_chunks = ones(size(transitions,1),2);
for i = 1:size(transitions,1)
    bypass_chunks(i,2) = transitions(i);
end
for i = 1:size(transitions,1)-1
    bypass_chunks(i+1,1) = transitions(i)+1;
end

% OUTLIER REMOVAL for Bypass Measurements
% Remove outliers more than three local scaled MAD from the local median
% Note: movmean (removing outliers more than three local standard deviations
% from the local mean) is too strongly affected by outliers

figure,plot(Bypass_datetime,Bypass_HCHO,'.','MarkerSize',15)
hold on
outlier_logical_array = [];
for i = 1:size(bypass_chunks,1)
    temp = isoutlier(Bypass_HCHO(bypass_chunks(i,1):bypass_chunks(i,2)),'gesd');
    outlier_logical_array = [outlier_logical_array; temp];
end
for j=1:length(outlier_logical_array)
    if outlier_logical_array(j)==1
        Bypass_HCHO(j) = NaN;
    end
end
plot(Bypass_datetime,Bypass_HCHO,'.','MarkerSize',15)
title('Bypass Outlier Removal')

% Remove NaNs
Bypass_posixtime = posixtime(Bypass_datetime);
temp_matrix = [Bypass_posixtime Bypass_HCHO];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Bypass_posixtime = temp_matrix(:,1);
Bypass_HCHO = temp_matrix(:,2);
Bypass_datetime = datetime(Bypass_posixtime,'ConvertFrom','posixtime');

% Interpolate Flag2 (representing equilibrated steps) onto Bypass_datetime
% 30 represents times when steps are equilibrated; defined prior to experiment
Bypass_Flag2 = interp1(p.plant_datetime,p.Flag2,Bypass_datetime);
Bypass_Equilibrated_Indices = find(Bypass_Flag2==30);

% Chunk equilibrated indices to break into their respective steps
Bypass_Equilibrated_Chunks = chunker(Bypass_Equilibrated_Indices);

% Show Flag 2 on Data
if s.engplot
    figure,plot(Bypass_datetime,Bypass_HCHO)
    hold on
    plot(Bypass_datetime,Bypass_Flag2,'.','MarkerSize',20)
    title('Bypass HCHO with Flag2')
end

Bypass_HCHO_avg = [];
Bypass_HCHO_avg_std = [];

for i = 1:size(Bypass_Equilibrated_Chunks,1)       
    Bypass_HCHO_avg(i)      = mean(Bypass_HCHO(Bypass_Equilibrated_Chunks(i,1):Bypass_Equilibrated_Chunks(i,2)));
    Bypass_HCHO_avg_std(i)  = std(Bypass_HCHO(Bypass_Equilibrated_Chunks(i,1):Bypass_Equilibrated_Chunks(i,2)))/sqrt(size(Bypass_HCHO(Bypass_Equilibrated_Chunks(i,1):Bypass_Equilibrated_Chunks(i,2)),1));
end

% For each bypass step with a given average and std dev, generate n number 
% of random numbers where n corresponds to the number of chamber
% measurements in that step

bypass_random_pts = []; % Randomly generated points for the bypass
for i = 1:size(Chamber_Equilibrated_Chunks,1)
   temp = size(Chamber_Equilibrated_Chunks(i,1):Chamber_Equilibrated_Chunks(i,2),2);
   r = normrnd(Bypass_HCHO_avg(i),Bypass_HCHO_avg_std(i),[1,temp])';
   bypass_random_pts = [bypass_random_pts; r];
end


%% Export MAT File for Import into R
% No NaNs should appear in the data!

chamber_HCHO = Chamber_HCHO_avg';
chamber_HCHO_std = Chamber_HCHO_avg_std';
bypass_HCHO = Bypass_HCHO_avg';
bypass_HCHO_std = Bypass_HCHO_avg_std';
save('chamber_blank_output_R_StepsAveraged.mat','chamber_HCHO','chamber_HCHO_std',...
    'bypass_HCHO','bypass_HCHO_std','chamber_HCHO_pts','chamber_groups','bypass_random_pts')


%% Commented Out Code

%     % Run Kolmogrov-Smirnov Tests to test for normality. The null hypothesis
%     % (Ho) is that the distribution is normal.
%     % h=0 accepts the null hypothesis (i.e., the data is normal). This corresponds to p > 0.05
% 
%     h_kstest_chamber = [];
%     h_kstest_bypass  = [];
%     p_kstest_chamber = [];
%     p_kstest_bypass  = [];
%     
%     for i=1:size(equilibrated_chunks,1)
%         kstest_chamber = Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2));
%         kstest_bypass = bypass_HCHO_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2));
%         [h_kstest_chamber(i),p_kstest_chamber(i)] = kstest(kstest_chamber);
%         [h_kstest_bypass(i),p_kstest_bypass(i)] = kstest(kstest_bypass);
%         figure,histogram(kstest_chamber,12)
%     end

% %On April 16, we decided not to perform t-tests on the chamber vs bypass
% %since removing points that were deemed statistically insignificant really
% %didn't change the fits and this procedure would then be hard to justify
% %in a future manuscript or talk
%   p_values = [];   
%     % Run t-tests to see if there's significant differences between bypass
%     % and chamber measurements. The null hypothesis (Ho) is that there's no
%     % significant difference between the chamber and bypass HCHO.
%     j = 1;
%     for i=1:size(equilibrated_chunks,1)
%         ttest_chamber = Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2));
%         ttest_bypass = bypass_HCHO_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2));
%         [h,pval] = ttest(ttest_chamber,ttest_bypass);
%         p_values(i) = pval;
%         
%         if h % h=1 rejects the null hypothesis (i.e., the bypass and chamber are different). This corresponds to p < 0.05
%             o.Chamber_HCHO(j) = mean(Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
%             o.Bypass_HCHO(j) = mean(bypass_HCHO_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
%             o.Chamber_HCHO_std(j) = std(Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
%             o.Bypass_HCHO_std(j) = std(bypass_HCHO_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
%             j = j+1;
%         end
%     end

% 90 s
% chamber_HCHO = o.Chamber_HCHO_90s;
% chamber_HCHO_std(1:size(chamber_HCHO,1)) = 0.018;
% bypass_HCHO = o.Bypass_HCHO_90s;
% bypass_HCHO_std(1:size(chamber_HCHO,1)) = 0.018;
% save('chamber_blank_output_R_SingleStep.mat','chamber_HCHO','chamber_HCHO_std',...
%     'bypass_HCHO','bypass_HCHO_std')

% o.Chamber_HCHO_90s = Chamber_HCHO(equilibrated_chunks(5,1):equilibrated_chunks(5,2));
% o.Bypass_HCHO_90s = bypass_HCHO_interp(equilibrated_chunks(5,1):equilibrated_chunks(5,2));
