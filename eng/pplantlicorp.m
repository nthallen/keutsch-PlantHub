function pplantlicorp(varargin);
% pplantlicorp( [...] );
% LICOR Pressure
h = timeplot({'LICOR_P'}, ...
      'LICOR Pressure', ...
      'Pressure', ...
      {'LICOR\_P'}, ...
      varargin{:} );
