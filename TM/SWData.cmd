%{
  #include "SWData.h"
  #ifdef SERVER
    SWData_t SWData;
  #endif
%}

%INTERFACE <SWData:DG/data>

&command
  : &SWTM * { if_SWData.Turf(); }
  ;
&SWTM
  : SWStat &SWStat { SWData.SWStat = $2; }
  : Flag &Flag { SWData.Flag = $2; }
  ;
&SWStat <unsigned char>
  : CO2 Control Idle { $0 = SWS_PCTRL_IDLE; }
  : Set %d { $0 = $2; }
  : CO2 Control Activate { $0 = SWS_PCTRL_ACTIVE; }
  : Time Warp { $0 = SWS_TIMEWARP; }
  : Shutdown { $0 = SWS_SHUTDOWN; }
  ;
&Flag <unsigned char>
  : Set %d { $0 = $2; }
  ;
