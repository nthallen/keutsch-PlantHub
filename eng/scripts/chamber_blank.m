%% Preamble

% Use this script for general workup of plant experiment data

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

%% Calculate VPD (Vapor Pressure Deficit) for Future Use

% From the twin cuvette paper, use 'crotch' equation (Goff-Gratch Equation)

% Saturation Water Vapor Pressure (hPa)
Tst = 373.15; %K Steam-point temperature
est = 1013.25; %hPa Steam-point pressure
T = p.TS_Temp + 273.15;%K
sat_H2O_vapor = zeros(size(T));

for i=1:length(T)
sat_H2O_vapor(i) = 10.^(-7.90298*((Tst/T(i))-1) + 5.02808*log10(Tst/T(i))...
    -1.3816e-7*(10.^(11.344*(1-(T(i)/Tst)))-1)...
    + 8.1328e-3*(10.^(-3.49149*((Tst/T(i))-1))-1)+log10(est));
end

% Calculate Vapor Pressure Deficit
P_chamber = 101.325; %kPa (must be in this unit)
p.VPD = 100*sat_H2O_vapor./(P_chamber*1000) - p.H2O_ppth*(10^-3);
% Will parse times corresponding to chamber sampling in the next section.


%% Parsing of Bypass vs Chamber HCHO
% Flag 51 - Sampling Bypass
% Flag 50 - Sampling Chamber
% Flag3 71 - Sampling LICOR
% Flag3 70 - Sampling PTR3

BypassIndices  = find(p.Flag1==51);
ChamberIndices = find(p.Flag1==50);
LICORIndices   = find(p.Flag3==71);
FILIFIndices    = find(p.Flag3==70);

% Parse FILIF Data
bypass_indices_FILIF = intersect(BypassIndices,FILIFIndices);
[bypass_indices_FILIF, ~] = RemoveFirstPoints(bypass_indices_FILIF, 10);
b.datetime_FILIF   = p.plant_datetime(bypass_indices_FILIF);
b.posixtime_FILIF  = p.Tplanteng_1(bypass_indices_FILIF);
b.hcho             = p.hcho(bypass_indices_FILIF);

chamber_indices_FILIF = intersect(ChamberIndices,FILIFIndices);
[chamber_indices_FILIF, ~] = RemoveFirstPoints(chamber_indices_FILIF, 180);
c.datetime_FILIF   = p.plant_datetime(chamber_indices_FILIF);
c.posixtime_FILIF  = p.Tplanteng_1(chamber_indices_FILIF);
c.hcho             = p.hcho(chamber_indices_FILIF);


% Separate bypass data from chamber data. All times in p are at 1 Hz.

% Bypass_datetime_1Hz  = p.plant_datetime(BypassIndices);
% Bypass_posixtime_1Hz = p.Tplanteng_1(BypassIndices);
% Bypass_CO2_1Hz       = p.CO2_ppm(BypassIndices);
% Bypass_H2O_1Hz       = p.H2O_ppth(BypassIndices);
% Bypass_HCHO_1Hz      = p.hcho(BypassIndices);
% 
% Chamber_datetime_1Hz  = p.plant_datetime(ChamberIndices);
% Chamber_posixtime_1Hz = p.Tplanteng_1(ChamberIndices);
% Chamber_CO2_1Hz       = p.CO2_ppm(ChamberIndices);
% Chamber_H2O_1Hz       = p.H2O_ppth(ChamberIndices);
% Chamber_HCHO_1Hz      = p.hcho(ChamberIndices);
% Chamber_VPD_1Hz       = p.VPD(ChamberIndices);

% Averaging
time_avg = 10; % Averaging time for data in seconds

if s.engplot
    figure,plot(b.datetime_FILIF,b.hcho,'.','MarkerSize',20)
    hold on
    plot(c.datetime_FILIF,c.hcho,'.','MarkerSize',20)
end

%% CHAMBER HCHO

[Chamber_posixtime, Chamber_HCHO] = binavg_plant(c.posixtime_FILIF,c.hcho,time_avg);

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
%Chamber_posixtime = posixtime(Chamber_datetime);
%temp_matrix = [Chamber_posixtime Chamber_HCHO];
%temp_matrix(any(isnan(temp_matrix),2),:) = [];
%Chamber_posixtime = temp_matrix(:,1);
%Chamber_HCHO = temp_matrix(:,2);
%Chamber_datetime = datetime(Chamber_posixtime,'ConvertFrom','posixtime');

Chamber_HCHO_avg = [];
Chamber_HCHO_avg_std = [];

for i = 1:size(chamber_chunks,1)       
    Chamber_HCHO_avg(i)      = nanmean(Chamber_HCHO(chamber_chunks(i,1):chamber_chunks(i,2)));
    Chamber_HCHO_avg_std(i)  = nanstd(Chamber_HCHO(chamber_chunks(i,1):chamber_chunks(i,2)))/sqrt(size(Chamber_HCHO(chamber_chunks(i,1):chamber_chunks(i,2)),1));
end

%% BYPASS HCHO
% First, obtain averaged bypass HCHO for each step of an experimental run with
% corresponding standard deviation of the mean

[Bypass_posixtime, Bypass_HCHO] = binavg_plant(b.posixtime_FILIF,b.hcho,time_avg);

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


% Chunk the Outlier-Removed Bypass Data to Prepare for Averaging
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

Bypass_HCHO_avg = [];
Bypass_HCHO_avg_std = [];

for i = 1:size(bypass_chunks,1)       
    Bypass_HCHO_avg(i)      = mean(Bypass_HCHO(bypass_chunks(i,1):bypass_chunks(i,2)));
    Bypass_HCHO_avg_std(i)  = std(Bypass_HCHO(bypass_chunks(i,1):bypass_chunks(i,2)))/sqrt(size(Bypass_HCHO(bypass_chunks(i,1):bypass_chunks(i,2)),1));
end

%% Export MAT File for Import into R
% No NaNs should appear in the data!

chamber_HCHO = Chamber_HCHO_avg';
chamber_HCHO_std = Chamber_HCHO_avg_std';
bypass_HCHO = Bypass_HCHO_avg';
bypass_HCHO_std = Bypass_HCHO_avg_std';
save('chamber_blank_output_R_StepsAveraged.mat','chamber_HCHO','chamber_HCHO_std',...
    'bypass_HCHO','bypass_HCHO_std')
