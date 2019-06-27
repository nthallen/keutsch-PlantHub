%% Insertion or Modification of Flag2

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
Plant_dir = ['D:\Plant\RAW\',s.plant_run_date,'\'];
addpath(Plant_dir)

% Load data files
p = load('planteng_1.mat'); %Loads all variables from plant computer to the 'p' structure

%% Time

% Convert Plant Posixtime to Datetime
p.plant_datetime = datetime(p.Tplanteng_1,'ConvertFrom','posixtime');

% Convert to local time if specified in config.ini file
if s.local_time_convert
    p.plant_datetime = p.plant_datetime - hours(s.time_adjust);
end

%% Define Times for Flag 2
% Set Flag2 = 30

p.Flag2 = zeros(size(p.plant_datetime,1),1);

% Step 1
start = find(p.plant_datetime.Hour==9 & p.plant_datetime.Minute==45);
final = find(p.plant_datetime.Hour==11 & p.plant_datetime.Minute==00);
p.Flag2(start(1):final(1)) = 30;

% Step 2
start = find(p.plant_datetime.Hour==12 & p.plant_datetime.Minute==5);
final = find(p.plant_datetime.Hour==13 & p.plant_datetime.Minute==20);
p.Flag2(start(1):final(1)) = 30;

% Step 3
start = find(p.plant_datetime.Hour==14 & p.plant_datetime.Minute==25);
final = find(p.plant_datetime.Hour==15 & p.plant_datetime.Minute==40);
p.Flag2(start(1):final(1)) = 30;

% Step 4
start = find(p.plant_datetime.Hour==16 & p.plant_datetime.Minute==45);
final = find(p.plant_datetime.Hour==18 & p.plant_datetime.Minute==00);
p.Flag2(start(1):final(1)) = 30;

% Step 5
start = find(p.plant_datetime.Hour==19 & p.plant_datetime.Minute==5);
final = find(p.plant_datetime.Hour==20 & p.plant_datetime.Minute==20);
p.Flag2(start(1):final(1)) = 30;

% Step 6
start = find(p.plant_datetime.Hour==21 & p.plant_datetime.Minute==25);
final = find(p.plant_datetime.Hour==22 & p.plant_datetime.Minute==40);
p.Flag2(start(1):final(1)) = 30;

%% Resave planteng_1.mat
clear 'p.plant_datetime'
save('planteng_1.mat','-struct','p')