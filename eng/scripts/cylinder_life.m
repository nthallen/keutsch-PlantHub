%% Preamble
% Determine the amount of time a gas will last in a cylinder assuming a
% continuous flow being controlled by a mass flow controller (MFC).

% We're assuming that 298 K is the temperature at which the cylinder was
% filled by the gas company

% The volumes of common Airgas cylinders (in L) can be found at
% https://ehs.mit.edu/site/sites/default/files/files/Airgas%20Cylinder%20info%20charts.pdf
% Martha has a volume of 4.13 L and is filled to 50 psi for safety reasons

% Time is in days

%% Initialization

flow_rate = 40; % sccm
initial_cylinder_pressure = 2500; % psi
cylinder_volume = 49.8; % L
standard_pressure = 14.7; % psi

%% Lifetime

t = (initial_cylinder_pressure*cylinder_volume)/(standard_pressure*(flow_rate/1000))/(60*24);

fprintf('The lifetime (in days) of the cylinder assuming continuous gas flow is = %f\n', t); 
