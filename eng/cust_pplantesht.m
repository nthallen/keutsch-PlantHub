function cust_pplantesht(h)
% cust_pplantesht(h)
% Customize plot created by pplantesht

% pplantesht's definition:

% function pplantesht(varargin);
% % pplantesht( [...] );
% % Eng SHT state
% h = timeplot({'SHT31_state'}, ...
%       'Eng SHT state', ...
%       'SHT state', ...
%       {'SHT31\_state'}, ...
%       varargin{:} );

% Example customizations include:
set(h,'LineStyle','none','Marker','.');
%   ax = get(h(1),'parent');
%   set(ax,'ylim',[0 800]);
