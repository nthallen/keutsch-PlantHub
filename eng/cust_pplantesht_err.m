function cust_pplantesht_err(h)
% cust_pplantesht_err(h)
% Customize plot created by pplantesht_err

% pplantesht_err's definition:

% function pplantesht_err(varargin);
% % pplantesht_err( [...] );
% % Eng SHT err
% h = timeplot({'SHT31_I2C_err'}, ...
%       'Eng SHT err', ...
%       'SHT err', ...
%       {'SHT31\_I2C\_err'}, ...
%       varargin{:} );

% Example customizations include:
  set(h,'LineStyle','none','Marker','.');
%  ax = get(h(1),'parent');
%  set(ax,'ylim',[0 800]);
