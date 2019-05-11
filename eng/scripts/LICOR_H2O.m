% LICOR_H2O.m
%
%% Initialization

T_bub = 273.15 + 23.8; %K
P_bub = 982; %hPa
T_mix = 273.15 + 23.0; %K 296.25 redundant
P_mix = 995; %hPa redundant
WetFlow = 200; %sccm
DryFlow = 300; %sccm

%% Equations


% Determine Saturation Mixing Ratio over Liquid Water using
% Murphy & Koop relationsclehip (SMRwmk in ppmv)
% [SMRwmk,SMRimk]=smrmk(T of bubbler in K, P of bubbler in hPa);
[SMRwmk,SMRimk] = smrmk(T_bub,P_bub);

% Determine Number Density of Mixed Flow (M in #/cc)
% (Use with T of mixed flow in K, P of mixed flow in hPa)
M_mix = 6.02e23/8.314/1e6*100*P_mix/T_mix;

% Determine Number Density of Water Vapor in Mixed Flow (#/cc)
H2Oadd_Bub=M_mix.*(SMRwmk*1e-6).*(WetFlow./(WetFlow+DryFlow));

% Determine Mixing Ratio of Water Vapor in Mixed Flow (ppmv)
VMRadd_Bub=H2Oadd_Bub./M_mix*1e6;

disp({'The mixing ratio of H2O (ppmv) is ',num2str(VMRadd_Bub)} )