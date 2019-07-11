clear variables
close all

basePath = 'D:\Plant\RAW\190708.3\';



%% import HLX manager data:

path = strcat(basePath,'Manager\');

files = dir(path);
x = {files.name};
A.Ozone = [];
A.Time = [];
A.RH = [];
A.T = [];
for i=3:size(files,1)
    filename = strcat(path,x{1,i});
    info = h5info(filename);
    att = {info.Attributes.Value};
    fileTime = att{10}; % Matlab time local
    Gases = h5read(filename,'/Gases/data/');
    A.T = [A.T; Gases.Temperature0x2Dact];
    A.RH = [A.RH; Gases.RH0x2Dact];
    A.Ozone = [A.Ozone; Gases.Ozone0x2Dact];
    A.Time = [A.Time; Gases.relativeTime0x5Bs0x5D/24/3600 + fileTime];
end
MGR = A;
clear att filename files fileTime Gases i info x A

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

% We may have to play around with the primary ion below. 37 was an isotope

[~, index] = min(abs(HLX.MassList - 37.042));
m36 = HLX.Cps(:,index)/0.00401*sqrt(100/37.042); % This is an isotope abundance norm; duty cycle accounted for
[~, index] = min(abs(HLX.MassList - 54.055));
m54 = HLX.Cps(:,index)*sqrt(100/54.055);
HLX.primIons = m36+m54;
HLX.primIonsSm = smooth(HLX.primIons,25); % Smooth primary ion before normalizing

figure ('Color','white');
plot(HLX.time,HLX.primIons);
hold on
plot(HLX.time,HLX.primIonsSm);

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

% MVK+NH4 = 88.07624
% 

[~, MVKindex] = min(abs(HLX.MassList - 88.07624));
figure,plot(HLX.time,HLX.ndcps(:,MVKindex))