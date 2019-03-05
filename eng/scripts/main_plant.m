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

%% TIME AVERAGE HCHO DATA

[FILIF.time_1Hz, FILIF.hcho_1Hz] = binavg(FILIF.posixtime,FILIF.hcho,1);
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

%% Find Indices for each HCHO Step for Equilibration Correction

% 3600 s chosen since chamber equilibrates within this time
%p.hcho([Step1(1:3600); Step2(1:3600); Step3(1:3600)]) = NaN;

% Find Bypass Indices for each Step 
 Step1_bypass = find(p.FC0_Set == 0.22 & p.Flag == 51);
 Step2_bypass = find(p.FC0_Set == 0.44 & p.Flag == 51);
 Step3_bypass = find(p.FC0_Set == 0.89 & p.Flag == 51);
 Step4_bypass = find(p.FC0_Set == 1.32 & p.Flag == 51);
 Step5_bypass = find(p.FC0_Set == 1.7799 & p.Flag == 51);
 Step6_bypass = find(p.FC0_Set == 2.2199 & p.Flag == 51);

% Find Chamber Indices for each Step 
 Step1_chamber = find(p.FC0_Set == 0.22 & p.Flag == 50);
 Step2_chamber = find(p.FC0_Set == 0.44 & p.Flag == 50);
 Step3_chamber = find(p.FC0_Set == 0.89 & p.Flag == 50);
 Step4_chamber = find(p.FC0_Set == 1.32 & p.Flag == 50);
 Step5_chamber = find(p.FC0_Set == 1.7799 & p.Flag == 50);
 Step6_chamber = find(p.FC0_Set == 2.2199 & p.Flag == 50);
 
% Apply Bypass Indices to HCHO data and Extract Time for Steps 1 - 6
Step1_bypass_hcho = p.hcho(Step1_bypass);
Step1_bypass_time = posixtime(p.plant_datetime(Step1_bypass));
Step1_bypass_time = Step1_bypass_time - Step1_bypass_time(1); 

Step2_bypass_hcho = p.hcho(Step2_bypass);
Step2_bypass_time = posixtime(p.plant_datetime(Step2_bypass));
Step2_bypass_time = Step2_bypass_time - Step2_bypass_time(1);
 
Step3_bypass_hcho = p.hcho(Step3_bypass);
Step3_bypass_time = posixtime(p.plant_datetime(Step3_bypass));
Step3_bypass_time = Step3_bypass_time - Step3_bypass_time(1);

Step4_bypass_hcho = p.hcho(Step4_bypass);
Step4_bypass_time = posixtime(p.plant_datetime(Step4_bypass));
Step4_bypass_time = Step4_bypass_time - Step4_bypass_time(1); 

Step5_bypass_hcho = p.hcho(Step5_bypass);
Step5_bypass_time = posixtime(p.plant_datetime(Step5_bypass));
Step5_bypass_time = Step5_bypass_time - Step5_bypass_time(1);
 
Step6_bypass_hcho = p.hcho(Step6_bypass);
Step6_bypass_time = posixtime(p.plant_datetime(Step6_bypass));
Step6_bypass_time = Step6_bypass_time - Step6_bypass_time(1);



% Apply Chamber Indices to HCHO data and Extract Time for Steps 1 - 6 
Step1_chamber_hcho = p.hcho(Step1_chamber);
Step1_chamber_time = posixtime(p.plant_datetime(Step1_chamber));
Step1_chamber_time = Step1_chamber_time - Step1_chamber_time(1); 

Step2_chamber_hcho = p.hcho(Step2_chamber);
Step2_chamber_time = posixtime(p.plant_datetime(Step2_chamber));
Step2_chamber_time = Step2_chamber_time - Step2_chamber_time(1);
 
Step3_chamber_hcho = p.hcho(Step3_chamber);
Step3_chamber_time = posixtime(p.plant_datetime(Step3_chamber));
Step3_chamber_time = Step3_chamber_time - Step3_chamber_time(1);

Step4_chamber_hcho = p.hcho(Step4_chamber);
Step4_chamber_time = posixtime(p.plant_datetime(Step4_chamber));
Step4_chamber_time = Step4_chamber_time - Step4_chamber_time(1); 

Step5_chamber_hcho = p.hcho(Step5_chamber);
Step5_chamber_time = posixtime(p.plant_datetime(Step5_chamber));
Step5_chamber_time = Step5_chamber_time - Step5_chamber_time(1);
 
Step6_chamber_hcho = p.hcho(Step6_chamber);
Step6_chamber_time = posixtime(p.plant_datetime(Step6_chamber));
Step6_chamber_time = Step6_chamber_time - Step6_chamber_time(1);

% Save Data to a .mat file to send to IGOR for Equilibration Time Determination 
 save('HCHO_Steps.mat','Step1_bypass_time','Step1_bypass_hcho','Step1_chamber_time','Step1_chamber_hcho',...%,...
      'Step2_bypass_time','Step2_bypass_hcho','Step2_chamber_time','Step2_chamber_hcho',...
      'Step3_bypass_time','Step3_bypass_hcho','Step3_chamber_time','Step3_chamber_hcho',...
      'Step4_bypass_time','Step4_bypass_hcho','Step4_chamber_time','Step4_chamber_hcho',...
      'Step5_bypass_time','Step5_bypass_hcho','Step5_chamber_time','Step5_chamber_hcho',...
      'Step6_bypass_time','Step6_bypass_hcho','Step6_chamber_time','Step6_chamber_hcho');

% You may want to take these individual steps and now plot in Igor for
% equilibration times

%% Parsing of Bypass vs Chamber
% Flag 51 - Sampling Bypass
% Flag 50 - Sampling Chamber

BypassIndices = find(p.Flag==51);
ChamberIndices = find(p.Flag==50);

[BypassIndices, BypassIndices_rm]  = RemoveFirstPoints(BypassIndices, 60);  % MAY NEED TO CHANGE
[ChamberIndices, ChamberIndices_rm]  = RemoveFirstPoints(ChamberIndices, 90); % MAY NEED TO CHANGE

% Separate bypass data from chamber data

Bypass_datetime  = p.plant_datetime(BypassIndices);
Bypass_posixtime = p.Tplanteng_1(BypassIndices);
Bypass_CO2       = p.CO2_ppm(BypassIndices);
Bypass_H2O       = p.H2O_ppth(BypassIndices);
Bypass_HCHO      = p.hcho(BypassIndices);

Chamber_datetime  = p.plant_datetime(ChamberIndices);
Chamber_posixtime = p.Tplanteng_1(ChamberIndices);
Chamber_CO2       = p.CO2_ppm(ChamberIndices);
Chamber_H2O       = p.H2O_ppth(ChamberIndices);
Chamber_HCHO      = p.hcho(ChamberIndices);

if s.engplot
 figure
 subplot(2,2,1)
 plot(Bypass_datetime,Bypass_CO2)
 hold on
 plot(Chamber_datetime,Chamber_CO2)
 title('CO_2')
 subplot(2,2,2)
 plot(Bypass_datetime,Bypass_H2O)
 hold on
 plot(Chamber_datetime,Chamber_H2O)
 title('H_2O')
 subplot(2,2,3)
 plot(p.plant_datetime,p.SHT31_Temp)
 hold on
 plot(p.plant_datetime,p.TS_Temp)
 title('Temperature (Blue SHT, Orange Leaf)')% Leaf and Chamber Temp
 subplot(2,2,4)  % RH (Overlay with temp?)
 plot(p.plant_datetime,p.SHT31_RH)
 title('RH')
 
 [bypass_30s_time, bypass_30s_hcho] = binavg(Bypass_posixtime,Bypass_HCHO,30); % 
 [chamber_30s_time, chamber_30s_hcho] = binavg(Chamber_posixtime,Chamber_HCHO,30);
 bypass_30s_datetime = datetime(bypass_30s_time,'ConvertFrom','posixtime'); % Convert to datetime
 chamber_30s_datetime = datetime(chamber_30s_time,'ConvertFrom','posixtime');
 
 figure,plot(Bypass_datetime,Bypass_HCHO)
 hold on
 plot(Chamber_datetime,Chamber_HCHO)
 hold on
 plot(p.plant_datetime,5*p.FC0_Set)
 hold on
 plot(FILIF.datetime1Hz,FILIF.Flag)
end

%% Chunk and Interpolation of Bypass H2O Data
% Chunking bypass data in order to average each bypass chunk for later
% interpolation

bypass_chunks = chunker(BypassIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

l=size(bypass_chunks,1);
bypass_avg_chunks=nan(l,1);
time_avg_bypass=nan(1,1);

for i=1:l
   j = bypass_chunks(i,1):bypass_chunks(i,2);
   bypass_avg_chunks(i) = nanmean(p.H2O_ppth(j));
   time_avg_bypass(i) = nanmean(p.Tplanteng_1(j));
end

time_avg_bypass_datetime = datetime(time_avg_bypass,'ConvertFrom','posixtime');

% Interpolation of bypass data
chamber_bypass_H2O = interp1(time_avg_bypass_datetime,bypass_avg_chunks,Chamber_datetime);


figure,plot(Chamber_datetime,Chamber_H2O,'.')
hold on
plot(Chamber_datetime,chamber_bypass_H2O,'.')

% Obtain difference between chamber and bypass H2O

H2O_diff = Chamber_H2O - chamber_bypass_H2O;

figure,plot(Chamber_datetime,H2O_diff,'.')

%% Define P_in and P_out

%p.p_in = (p.hcho + s.source_term)*(1-s.wall_loss_fraction);
%p.p_out = p.hcho/(1-s.wall_loss_fraction) - s.source_term;

p.p_in = p.hcho*(1-s.wall_loss_fraction) + s.source_term;
p.p_out = p.hcho;

%% Chunk and Interpolation of Bypass and Chamber HCHO Data
% Chunking bypass and chamber data

p_in_chunks = chunker(BypassIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

l=size(p_in_chunks,1);
p_in_avg_chunks=nan(l,1);
time_p_in_avg_chunks=nan(1,1);

for i=1:l
   j = p_in_chunks(i,1):p_in_chunks(i,2);
   p_in_avg_chunks(i) = nanmean(p.p_in(j));
   time_p_in_avg_chunks(i) = nanmean(p.Tplanteng_1(j));
end

datetime_p_in_avg_chunks = datetime(time_p_in_avg_chunks,'ConvertFrom','posixtime');

figure,plot(datetime_p_in_avg_chunks,p_in_avg_chunks,'.')

%% Chunk and Interpolation of Bypass and Chamber HCHO Data
% Chunking bypass and chamber data

p_out_chunks = chunker(ChamberIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

l=size(p_out_chunks,1);
p_out_avg_chunks=nan(l,1);
time_p_out_avg_chunks=nan(1,1);

for i=1:l
   j = p_out_chunks(i,1):p_out_chunks(i,2);
   p_out_avg_chunks(i) = nanmean(p.p_out(j));
   time_p_out_avg_chunks(i) = nanmean(p.Tplanteng_1(j));
end

datetime_p_out_avg_chunks = datetime(time_p_out_avg_chunks,'ConvertFrom','posixtime');

figure
plot(datetime_p_in_avg_chunks,p_in_avg_chunks,'.-')
hold on
plot(datetime_p_out_avg_chunks,p_out_avg_chunks,'.-')

%% Calculation of Flux

% Interpolation of p_in data. p_in_interpol is used for flux
p.p_in_interpol = interp1(datetime_p_in_avg_chunks,p_in_avg_chunks,p.plant_datetime);

figure,plot(p.plant_datetime,p.p_in_interpol,'.')

% P_out - P_in Difference
p.diff = p.p_out - p.p_in_interpol;
f.diff = p.diff(ChamberIndices); % Still in units of ppbv
f.diff_datetime = p.plant_datetime(ChamberIndices);
f.posixtime = posixtime(f.diff_datetime);

% Calculation of Number Density in Air
s.pressure = 1; %atm
s.temperature = nanmean(p.SHT31_Temp) + 273.15; %K
s.R = 82.06; % cm^3*atm / mol*K
n_air = s.pressure*(6.02*10^23)/(s.R*s.temperature);

% Convert HCHO from ppbv to molecules/cm^3
f.n_hcho = f.diff * n_air * (10^-9);

f.flux = (s.flow_rate/s.leaf_area)*f.n_hcho; % molecules/cm^2 * min

f.flux_final = f.flux*(100)^2 *(6.02*10^-23)*(1/60)*(10^9); % nmol/m^2*s


[f.flux_avg_time, f.flux_avg_hcho] = binavg(f.posixtime,f.flux_final,120);
f.flux_avg_datetime = datetime(f.flux_avg_time,'ConvertFrom','posixtime');

figure,plot(f.flux_avg_datetime,f.flux_avg_hcho)
