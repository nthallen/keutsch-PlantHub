function cust_pplantets_error(h)
% cust_pplantets_error(h)
% Customize plot created by pplantets_error

% pplantets_error's definition:

% function pplantets_error(varargin);
% % pplantets_error( [...] );
% % Eng TS error
% h = timeplot({'TS_error'}, ...
%       'Eng TS error', ...
%       'TS error', ...
%       {'TS\_error'}, ...
%       varargin{:} );

% Example customizations include:
set(h,'LineStyle','none','Marker','.');
%   ax = get(h(1),'parent');
%   set(ax,'ylim',[0 800]);
