function pplantlicors(varargin);
% pplantlicors( [...] );
% LICOR Status
h = ne_dstat({
  'Fresh', 'LICOR_Status', 0; ...
	'Overflow', 'LICOR_Status', 1 }, 'Status', varargin{:} );
