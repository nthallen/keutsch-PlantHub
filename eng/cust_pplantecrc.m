function cust_pplantecrc(h)
% cust_pplantecrc(h)
% Customize plot created by pplantecrc

% pplantecrc's definition:

% function pplantecrc(varargin);
% % pplantecrc( [...] );
% % Eng CRC err
% h = timeplot({'SHT31_CRC_err'}, ...
%       'Eng CRC err', ...
%       'CRC err', ...
%       {'SHT31\_CRC\_err'}, ...
%       varargin{:} );

% Example customizations include:
set(h,'LineStyle','none','Marker','.');
%   ax = get(h(1),'parent');
%   set(ax,'ylim',[0 800]);
