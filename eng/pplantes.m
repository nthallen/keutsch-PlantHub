function pplantes(varargin);
% pplantes( [...] );
% Eng States
h = timeplot({'I2C_state','I2C_ts_state','I2C_sht_state'}, ...
      'Eng States', ...
      'States', ...
      {'I2C\_state','I2C\_ts\_state','I2C\_sht\_state'}, ...
      varargin{:} );
