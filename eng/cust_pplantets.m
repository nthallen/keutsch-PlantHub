function cust_pplantets(h)
% cust_pplantets(h)
% Customize plot created by pplantets

% pplantets's definition:

% function pplantets(varargin);
% % pplantets( [...] );
% % Eng TS state
% h = timeplot({'TS_state'}, ...
%       'Eng TS state', ...
%       'TS state', ...
%       {'TS\_state'}, ...
%       varargin{:} );

% Example customizations include:
set(h,'LineStyle','none','Marker','.');
%   ax = get(h(1),'parent');
%   set(ax,'ylim',[0 800]);
