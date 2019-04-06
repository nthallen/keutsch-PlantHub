% smrmk.m
% function [SMRwmk,SMRimk]=smrmk(T,P);
% T [K], P [mbar]
function [SMRwmk,SMRimk]=smrmk(T,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -> Vapor over Liquid as a function of ambient T & P
% Muphy and Koop, 2005 
% ** Integrated from Clasius-Clapyeron, pinned to Triple-point Data **
logPw=54.842763-6763.22./T-4.210.*log(T)+0.000367.*T+...
    tanh(0.0415.*(T-218.8)).*...
    (53.878-1331.22./T-9.44523.*log(T)+0.014025.*T);
Pw=exp(logPw);
SMRwmk=1e4.*Pw./P;

% Bolton
% (source: Holger Voemel's website -> Bolton, 1980)
% NOTE: This is the equation used in Harvard water vapor calibrations
Pw=6.112*exp(17.67.*(T-273.16)./(T-273.16+243.5));
SMRwbl=1e6.*Pw./P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -> Vapor over Ice as a function of ambient T & P
% Muphy and Koop, 2005 
% ** Integrated from Clasius-Clapyeron, pinned to Triple-point Data **
Pi=exp(9.550426 - 5723.265./T + 3.53068.*log(T) - 0.00728332.*T);
SMRimk=1e4.*Pi./P;

% Marti and Mauersberger 
% ** This equation based on direct measurments of pH2O down to ~170 K **
% (source: their paper)
% T in [K], P in [Pa]
log10Pi=-2663.5./T+12.537;
Pi=10.^(log10Pi); % =exp(2.303*(-2663.5./T+12.537));
SMRimm=1e4*Pi./P;

