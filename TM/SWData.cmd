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
  : Valve Switching Activate { $0 = SWS_VALVE_SWITCHING_ACTIVE; }
  : Set %d { $0 = $2; }
  : Valve Switching Idle { $0 = SWS_VALVE_SWITCHING_IDLE; }
  : Auto Start Flow Idle { $0 = SWS_AUTO_START_FLOW_IDLE; }
  : Auto Start Flow Activate { $0 = SWS_AUTO_START_FLOW_ACTIVE; }
  : Auto Start Flow Execute { $0 = SWS_AUTO_START_FLOW_EXECUTE; }
  : FCC Test Idle { $0 = SWS_FCC_TEST_IDLE; }
  : FCC Test Execute { $0 = SWS_FCC_TEST; }
  : Time Warp { $0 = SWS_TIMEWARP; }
  : Shutdown { $0 = SWS_SHUTDOWN; }
  ;
&Flag <unsigned char>
  : Set %d { $0 = $2; }
  ;
