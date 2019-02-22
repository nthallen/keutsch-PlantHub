function fig = gplante(varargin);
% gplante(...)
% Eng
ffig = ne_group(varargin,'Eng','pplantetsa','pplantetsc','pplantets','pplantets_error','pplanteovf','pplantesht','pplantesht_err','pplantecrc','pplantes','pplanteb','pplantestatus');
if nargout > 0 fig = ffig; end
