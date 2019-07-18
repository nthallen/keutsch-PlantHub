%% Plant Workup Script for Use with PTR3 Data
% July 2019

%% SETTINGS AND FOLDER INITIALIZATION
% Settings for a run are located in a user-defined config.ini file. Please
% see an example config.ini file for more details.

%clear variables
%close all

% Load configuration file containing settings
s = ini2struct('config.ini');
fields = fieldnames(s); % List field names in the settings structure
for i = 1:numel(fields) 
    s.(fields{i}) = str2double(s.(fields{i})); % Convert character strings to numerical values
end

s.plant_run_date = num2str(s.plant_run_date);
Plant_dir = ['D:\Plant\RAW\',s.plant_run_date,'\'];
addpath(Plant_dir)

% Load plant data file
p = load('planteng_1.mat'); %Loads all variables from plant computer to the 'p' structure

% Convert Plant Posixtime to Datetime
p.plant_datetime = datetime(p.Tplanteng_1,'ConvertFrom','posixtime');

% Convert to local time if specified in config.ini file
if s.local_time_convert
    p.plant_datetime = p.plant_datetime - hours(s.time_adjust);
end

% Correct plant posixtime array for possible change due to local time conversion
p.Tplanteng_1 = posixtime(p.plant_datetime);

% basePath for PTR3 data
basePath = Plant_dir;

%% import HLX manager data:

% path = strcat(basePath,'Manager\');
% 
% files = dir(path);
% x = {files.name};
% A.Ozone = [];
% A.Time = [];
% A.RH = [];
% A.T = [];
% for i=3:size(files,1)
%     filename = strcat(path,x{1,i});
%     info = h5info(filename);
%     att = {info.Attributes.Value};
%     fileTime = att{10}; % Matlab time local
%     Gases = h5read(filename,'/Gases/data/');
%     A.T = [A.T; Gases.Temperature0x2Dact];
%     A.RH = [A.RH; Gases.RH0x2Dact];
%     A.Ozone = [A.Ozone; Gases.Ozone0x2Dact];
%     A.Time = [A.Time; Gases.relativeTime0x5Bs0x5D/24/3600 + fileTime];
% end
% MGR = A;
% clear att filename files fileTime Gases i info x A

%% import HLX stuff:
path = strcat(basePath);
filename = strcat(path,'_StickResult.hdf5');
HLX.Cps = hdf5read(filename,'/Cps');
HLX.MassList = hdf5read(filename,'/MassList');
%HLX.time = hdf5read(filename,'/MatlabTimes')-4/24;
HLX.time = hdf5read(filename,'/UnixTimestamps');
tmp = hdf5read(filename,'/ElementNames');
HLX.Composition = double(hdf5read(filename,'/ElementalComposition'));
HLX.Elements = cell(size(tmp,1),1);
for i = 1:size(tmp,1)
    HLX.Elements(i) = {tmp(i).Data};
end
clear tmp i

% Convert the PTR3 data to same time zone as plant data
HLX.datetime = datetime(HLX.time,'ConvertFrom','posixtime') - hours(s.time_adjust);

% We may have to play around with the primary ion below. 37 was an isotope

[~, index] = min(abs(HLX.MassList - 37.042));
m36 = HLX.Cps(:,index)/0.00401*sqrt(100/37.042); % This is an isotope abundance norm; duty cycle accounted for
[~, index] = min(abs(HLX.MassList - 54.055));
m54 = HLX.Cps(:,index)*sqrt(100/54.055);
%[~, index] = min(abs(HLX.MassList - 18.03437));
%m18 = HLX.Cps(:,index)/1*sqrt(100/18.03437);
HLX.primIons = m36+m54;
HLX.primIonsSm = smooth(HLX.primIons,25); % Smooth primary ion before normalizing

figure ('Color','white');
plot(HLX.datetime,HLX.primIons);
hold on
plot(HLX.datetime,HLX.primIonsSm);

% smooth Cps:
for i=1:size(HLX.MassList,1)
    HLX.Cps(:,i) = smooth(HLX.Cps(:,i),30); %Smooth for 30 s
    i=i;
end

% create nCps:
HLX.dcps = HLX.Cps*nan;
HLX.ndcps = HLX.Cps*nan;
for i=1:size(HLX.MassList,1)
    HLX.dcps(:,i) = HLX.Cps(:,i)*sqrt(100/HLX.MassList(i)); %Duty-cycle corrected cps
    HLX.ndcps(:,i) = HLX.dcps(:,i)*1e6./HLX.primIonsSm; %Normalize to 1 million primary ions
    i=i;
end

HLX.Cps(HLX.Cps<0) = nan;

%% Pull Out Masses of Interest

% Use SIS Isotopes website to find masses!

% MVK+NH4 or MACR+NH4  = 88.07624
% ISOPOOH+NH4          = 136.09738
% MEK+NH4              = 90.09189
% ISOPRENE+NH4         = 86.09697
% MethylSalicylate+NH4 = 170.08173
% Acetone+NH4          = 76.07624
% 2-methyl-2-vinyloxirane+NH4 = 102.09189

[~, m.MVKMACRindex] = min(abs(HLX.MassList - 88.07624));
m.MVKMACR = HLX.ndcps(:,m.MVKMACRindex);

[~, m.Acetoneindex] = min(abs(HLX.MassList - 76.07624));
m.Acetone = HLX.ndcps(:,m.Acetoneindex);

[~, m.ISOPOOHindex] = min(abs(HLX.MassList - 136.09738));
m.ISOPOOH = HLX.ndcps(:,m.ISOPOOHindex);

% Note: Of course, there could be other chemical structures for some of
% these masses, though a chemical name is listed to help with discussion.


%% Time Sync Plant and PTR3 Datasets
% Treat Plant datetime as the master clock and the PTR3 datetime is either
% ahead or behind this time (in seconds)

% Convert the PTR3 data to same time zone as plant data. This was copied in
% the event that someone accidentally runs this section multiple times
HLX.datetime = datetime(HLX.time,'ConvertFrom','posixtime') - hours(s.time_adjust);

file_exist_check = exist(fullfile(Plant_dir,'TimeOffset.mat'), 'file');

if file_exist_check == 2
    load('TimeOffset.mat')
else
    figure,plot(HLX.datetime,m.Acetone)
    hold on
    plot(p.plant_datetime,500*p.Flag3)
    
    prompt = 'What is the value of time_offset (in seconds)? ';
    time_offset = input(prompt);
end
 
if exist('time_offset','var') == 1
    save(fullfile(Plant_dir,'TimeOffset.mat'),'time_offset');
else
    disp('You have not made a time_offset variable (in seconds) to sync PTR3 and Plant datetimes!')
end

HLX.datetime = HLX.datetime + seconds(time_offset);

figure,plot(HLX.datetime,m.Acetone)
hold on
plot(p.plant_datetime,500*p.Flag3)

%% Move PTR3 Data to Plant Time Basis

% Now that the plant and PTR3 datetimes are synced, let's get the PTR3 on
% the plant datetime time basis for easier processing. We'll do this by
% interpolation

p.PTR3.Acetone = interp1(HLX.datetime,m.Acetone,p.plant_datetime);
p.PTR3.MVKMACR = interp1(HLX.datetime,m.MVKMACR,p.plant_datetime);
p.PTR3.ISOPOOH = interp1(HLX.datetime,m.ISOPOOH,p.plant_datetime);

%% Parsing of Bypass vs Chamber and LICOR vs PTR3 Sampling Periods
% Flag 51 - Sampling Bypass
% Flag 50 - Sampling Chamber
% Flag3 71 - Sampling LICOR
% Flag3 70 - Sampling PTR3

BypassIndices  = find(p.Flag==51);
ChamberIndices = find(p.Flag==50);
LICORIndices   = find(p.Flag3==71);
PTR3Indices    = find(p.Flag3==70);

% LICOR Data
bypass_indices_LICOR = intersect(BypassIndices,LICORIndices);
[bypass_indices_LICOR, ~] = RemoveFirstPoints(bypass_indices_LICOR, 35);
b.datetime_LICOR  = p.plant_datetime(bypass_indices_LICOR);
b.posixtime_LICOR = p.Tplanteng_1(bypass_indices_LICOR);
b.CO2_ppm         = p.CO2_ppm(bypass_indices_LICOR);
b.H2O_ppth        = p.H2O_ppth(bypass_indices_LICOR);

chamber_indices_LICOR = intersect(ChamberIndices,LICORIndices);
[chamber_indices_LICOR, ~] = RemoveFirstPoints(chamber_indices_LICOR, 35);
c.datetime_LICOR  = p.plant_datetime(chamber_indices_LICOR);
c.posixtime_LICOR = p.Tplanteng_1(chamber_indices_LICOR);
c.CO2_ppm         = p.CO2_ppm(chamber_indices_LICOR);
c.H2O_ppth        = p.H2O_ppth(chamber_indices_LICOR);

% PTR3 Data
bypass_indices_PTR3 = intersect(BypassIndices,PTR3Indices);
[bypass_indices_PTR3, ~] = RemoveFirstPoints(bypass_indices_PTR3, 35);
b.datetime_PTR3   = p.plant_datetime(bypass_indices_PTR3);
b.posixtime_PTR3  = p.Tplanteng_1(bypass_indices_PTR3);
b.Acetone         = p.PTR3.Acetone(bypass_indices_PTR3);
b.MVKMACR         = p.PTR3.MVKMACR(bypass_indices_PTR3);

chamber_indices_PTR3 = intersect(ChamberIndices,PTR3Indices);
[chamber_indices_PTR3, ~] = RemoveFirstPoints(chamber_indices_PTR3, 35);
c.datetime_PTR3   = p.plant_datetime(chamber_indices_PTR3);
c.posixtime_PTR3  = p.Tplanteng_1(chamber_indices_PTR3);
c.Acetone         = p.PTR3.Acetone(chamber_indices_PTR3);
c.MVKMACR         = p.PTR3.MVKMACR(chamber_indices_PTR3);


if s.engplot
    figure,plot(b.datetime_LICOR,b.H2O_ppth,'.','MarkerSize',20)
    hold on
    plot(c.datetime_LICOR,c.H2O_ppth,'.','MarkerSize',20)
end
