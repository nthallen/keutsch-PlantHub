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


%% Parsing of Bypass vs Chamber
% Flag 51 - Sampling Bypass
% Flag 50 - Sampling Chamber

BypassIndices = find(p.Flag==51);
ChamberIndices = find(p.Flag==50);

[BypassIndices, BypassIndices_rm]  = RemoveFirstPoints(BypassIndices, 60);
[ChamberIndices, ChamberIndices_rm]  = RemoveFirstPoints(ChamberIndices, 1000);

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
Chamber_VPD_1Hz       = p.VPD(ChamberIndices);

% Averaging
time_avg = 10; % Averaging time for data in seconds

[~, Bypass_CO2] = binavg_plant(Bypass_posixtime_1Hz, Bypass_CO2_1Hz,time_avg);
[~, Bypass_H2O] = binavg_plant(Bypass_posixtime_1Hz, Bypass_H2O_1Hz,time_avg);
[Bypass_posixtime, Bypass_HCHO] = binavg_plant(Bypass_posixtime_1Hz, Bypass_HCHO_1Hz,time_avg);

% Remove NaNs
temp_matrix = [Bypass_posixtime' Bypass_CO2 Bypass_H2O Bypass_HCHO];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Bypass_posixtime = temp_matrix(:,1);
Bypass_CO2 = temp_matrix(:,2);
Bypass_H2O = temp_matrix(:,3);
Bypass_HCHO = temp_matrix(:,4);
Bypass_datetime = datetime(Bypass_posixtime,'ConvertFrom','posixtime');

[~, Chamber_CO2] = binavg_plant(Chamber_posixtime_1Hz, Chamber_CO2_1Hz,time_avg);
[~, Chamber_H2O] = binavg_plant(Chamber_posixtime_1Hz, Chamber_H2O_1Hz,time_avg);
[~, Chamber_HCHO] = binavg_plant(Chamber_posixtime_1Hz, Chamber_HCHO_1Hz,time_avg);
[Chamber_posixtime, Chamber_VPD] = binavg_plant(Chamber_posixtime_1Hz, Chamber_VPD_1Hz,time_avg);

% Remove NaNs
temp_matrix = [Chamber_posixtime' Chamber_CO2 Chamber_H2O Chamber_HCHO Chamber_VPD];
temp_matrix(any(isnan(temp_matrix),2),:) = [];
Chamber_posixtime = temp_matrix(:,1);
Chamber_CO2 = temp_matrix(:,2);
Chamber_H2O = temp_matrix(:,3);
Chamber_HCHO = temp_matrix(:,4);
Chamber_VPD = temp_matrix(:,5);
Chamber_datetime = datetime(Chamber_posixtime,'ConvertFrom','posixtime');


if s.engplot
 figure
 subplot(2,2,1)
 plot(Bypass_datetime,Bypass_CO2,'.','MarkerSize',15)
 hold on
 plot(Chamber_datetime,Chamber_CO2,'.','MarkerSize',15)
 title('CO_2')
 subplot(2,2,2)
 plot(Bypass_datetime,Bypass_H2O,'.','MarkerSize',15)
 hold on
 plot(Chamber_datetime,Chamber_H2O,'.','MarkerSize',15)
 title('H_2O')
 subplot(2,2,3)
 plot(p.plant_datetime,p.SHT31_Temp)
 hold on
 plot(p.plant_datetime,p.TS_Temp)
 title('Temperature (Blue SHT, Orange Leaf)')% Leaf and Chamber Temp
 subplot(2,2,4)  % RH (Overlay with temp?)
 plot(p.plant_datetime,p.SHT31_RH)
 title('RH')
 

 figure,plot(Bypass_datetime,Bypass_HCHO,'.','MarkerSize',20)
 hold on
 plot(Chamber_datetime,Chamber_HCHO,'.','MarkerSize',20)
 hold on
 plot(p.plant_datetime,5*p.FC0_Set)
end

















%% Calculation of CO2 and H2O Differences Between Chamber and Bypass

% H2O Difference Calculation
% Chunk and Interpolation of Bypass H2O Data
bypass_chunks = chunker(BypassIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

len=size(bypass_chunks,1);
bypass_avg_chunks=nan(len,1);
time_avg_bypass=nan(len,1);

for i=1:len
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

f.H2O_diff = Chamber_H2O - chamber_bypass_H2O;

% CO2 Difference Calculation
% Chunk and Interpolation of Bypass CO2 Data
bypass_chunks = chunker(BypassIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

len=size(bypass_chunks,1);
bypass_avg_chunks=nan(len,1);
time_avg_bypass=nan(len,1);

for i=1:len
   j = bypass_chunks(i,1):bypass_chunks(i,2);
   bypass_avg_chunks(i) = nanmean(p.CO2_ppm(j));
   time_avg_bypass(i) = nanmean(p.Tplanteng_1(j));
end

time_avg_bypass_datetime = datetime(time_avg_bypass,'ConvertFrom','posixtime');

% Interpolation of bypass data
chamber_bypass_CO2 = interp1(time_avg_bypass_datetime,bypass_avg_chunks,Chamber_datetime);


figure,plot(Chamber_datetime,Chamber_CO2,'.')
hold on
plot(Chamber_datetime,chamber_bypass_CO2,'.')

% Obtain difference between chamber and bypass CO2

f.CO2_diff = Chamber_CO2 - chamber_bypass_CO2;

%% Calculate CO2 and H2O Fluxes

% Calculation of Molar Flow Rate (mol/s)
s.pressure = 1; %atm
s.temperature = 273.15; %K - Under standard conditions
s.R = 0.0821; % L*atm / mol*K
n_out = s.pressure*(s.flow_rate/1000)*(1/60)/(s.R*s.temperature); % mol/s

% Obtain CO2 and H2O Fluxes
f.flux_CO2 = f.CO2_diff*(10^-6)*n_out*(1/s.leaf_area)*(10^6)*(100)^2; % umol/m^2*s 
f.flux_H2O = f.H2O_diff*(10^-3)*n_out*(1/s.leaf_area)*(10^3)*(100)^2; % mmol/m^2*s

figure
plot(Chamber_datetime,f.flux_H2O)
title('H2O Flux')

figure
plot(Chamber_datetime,f.flux_CO2)
title('CO2 Flux')


%% Calculation of Stomatal Conductance

f.gs = f.flux_H2O./Chamber_VPD; % Bulk Stomatal Conductance

figure,plot(Chamber_datetime,f.gs)
title('H2O Stomatal Conductance')

%% Chunk and Interpolation of Bypass HCHO Data

bypass_HCHO_chunks = chunker(BypassIndices);

%Find average counts for each offline data chunk. Using interpolation 
%between the two adjacent offline counts averages, put onlines and average
%offlines on same time basis.

l=size(bypass_HCHO_chunks,1);
bypass_avg_HCHO_chunks=nan(l,1);
bypass_avg_time_chunks=nan(1,1);

for i=1:l
   j = bypass_HCHO_chunks(i,1):bypass_HCHO_chunks(i,2);
   bypass_avg_HCHO_chunks(i) = nanmean(p.hcho(j));
   bypass_avg_time_chunks(i) = nanmean(p.Tplanteng_1(j));
end

bypass_avg_datetime_chunks = datetime(bypass_avg_time_chunks,'ConvertFrom','posixtime');

figure,plot(bypass_avg_datetime_chunks,bypass_avg_HCHO_chunks,'.','MarkerSize',15)

%% Convert to Blank Chamber Out Using Bypass HCHO

c_out_blank = s.blank_conversion_slope*bypass_avg_HCHO_chunks + s.blank_conversion_intercept;

c_out_blank_interp = interp1(bypass_avg_datetime_chunks,c_out_blank,Chamber_datetime);

figure,plot(Chamber_datetime,c_out_blank_interp,'.','MarkerSize',15)

%% Calculation of HCHO Flux and Normalized HCHO Flux

% Take HCHO difference between c_out and c_out_blank
f.HCHO_diff = Chamber_HCHO - c_out_blank_interp;
f.flux_HCHO = f.HCHO_diff*(10^-9)*n_out*(1/s.leaf_area)*(10^9)*(100)^2; % nmol/m^2*s 

f.flux_HCHO_norm_gs = f.flux_HCHO./f.gs;

figure,plot(Chamber_datetime,f.flux_HCHO,'.')
hold on
plot(Chamber_datetime,f.flux_HCHO_norm_gs,'.')
hold on
plot(p.plant_datetime,p.Flag2)

figure,plot(Chamber_datetime,c_out_blank_interp,'.','MarkerSize',15)
hold on
plot(Chamber_datetime,Chamber_HCHO,'.','MarkerSize',15)
hold on
plot(p.plant_datetime,p.Flag2)


if isfield(p,'Flag2')
    
    Chamber_Flag2 = interp1(p.plant_datetime,p.Flag2,Chamber_datetime);
    Equilibrated_Indices  = find(Chamber_Flag2==30); % 30 represents times when steps are equilibrated; defined prior to experiment

    % Chunk equilibrated indices to break into their respective steps
    equilibrated_chunks = chunker(Equilibrated_Indices);
    
    o.flux_HCHO = [];
    o.flux_HCHO_norm_gs = [];
    o.Chamber_HCHO = [];
    o.flux_HCHO_std = [];
    o.flux_HCHO_norm_gs_std = [];
    o.Chamber_HCHO_std = [];
    o.Chamber_Blank_HCHO = [];
    o.Chamber_Blank_HCHO_std = [];
    
    for i = 1:size(equilibrated_chunks,1)
        o.flux_HCHO(i) = mean(f.flux_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.flux_HCHO_norm_gs(i) = mean(f.flux_HCHO_norm_gs(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.Chamber_HCHO(i) = mean(Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.Chamber_Blank_HCHO(i) = mean(c_out_blank_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        
        o.flux_HCHO_std(i) = std(f.flux_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.flux_HCHO_norm_gs_std(i) = std(f.flux_HCHO_norm_gs(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.Chamber_HCHO_std(i) = std(Chamber_HCHO(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
        o.Chamber_Blank_HCHO_std(i) = std(c_out_blank_interp(equilibrated_chunks(i,1):equilibrated_chunks(i,2)));
    end
end

%%

flux_HCHO = o.flux_HCHO';
flux_HCHO_norm_gs = o.flux_HCHO_norm_gs';
chamber_HCHO = o.Chamber_HCHO';
flux_HCHO_std = o.flux_HCHO_std';
flux_HCHO_norm_gs_std = o.flux_HCHO_norm_gs_std';
chamber_HCHO_std = o.Chamber_HCHO_std';
chamber_blank_HCHO = o.Chamber_Blank_HCHO';
chamber_blank_HCHO_std = o.Chamber_Blank_HCHO_std';
save('plant_output_R.mat','flux_HCHO','flux_HCHO_norm_gs','chamber_HCHO',...
    'flux_HCHO_std','flux_HCHO_norm_gs_std','chamber_HCHO_std',...
    'chamber_blank_HCHO','chamber_blank_HCHO_std')
save('plant_output_MATLAB.mat','o')

figure,plot(chamber_HCHO,flux_HCHO_norm_gs,'.','MarkerSize',15)

figure,plot(chamber_blank_HCHO,chamber_HCHO,'.','MarkerSize',15)
hold on
syms x
fplot(x,[0 20],'r','LineWidth',2)


%%

figure ('DefaultAxesFontSize',18)
% errorbar(Day1.chamber,Day1.flux,Day1.flux_err,Day1.flux_err,Day1.chamber_err,Day1.chamber_err,'b.','MarkerSize',12)
% hold on
% errorbar(Day2.chamber,Day2.flux,Day2.flux_err,Day2.flux_err,Day2.chamber_err,Day2.chamber_err,'r.','MarkerSize',12)
% hold on
% errorbar(Day3.chamber,Day3.flux,Day3.flux_err,Day3.flux_err,Day3.chamber_err,Day3.chamber_err,'g.','MarkerSize',12)

errorbar(chamber,flux,flux_err,flux_err,chamber_err,chamber_err,'.','MarkerSize',15, 'LineWidth',1)
hold on
syms x
fplot(-0.000537399976615926*x+0.000102510481285949,[0 20],'r','LineWidth',2)

title('Compensation Points for Plant F Golden Pothos')
xlabel('Chamber HCHO / ppbv')
ylabel('Normalized Flux')

%%
figure, errorbar(Day1.bypass,Day1.chamber,Day1.chamber_err,Day1.chamber_err,Day1.bypass_err,Day1.bypass_err,'.','MarkerSize',15)
hold on
errorbar(Day2.bypass,Day2.chamber,Day2.chamber_err,Day2.chamber_err,Day2.bypass_err,Day2.bypass_err,'.','MarkerSize',15)
syms x
fplot(0.961*x +0.25,[0 20],'r','LineWidth',2)
