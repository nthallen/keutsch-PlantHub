function fig = gplanttm(varargin);
% gplanttm(...)
% T Mbase
ffig = ne_group(varargin,'T Mbase','pplanttmtd','pplanttmcpu','pplanttmram','pplanttmd');
if nargout > 0 fig = ffig; end
