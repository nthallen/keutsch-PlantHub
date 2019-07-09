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
  : Flag2 &Flag2 { SWData.Flag2 = $2; }
  ;
&SWStat <unsigned char>
  : Valve Switching Activate { $0 = SWS_VALVE_SWITCHING_ACTIVE; }
  : Set %d { $0 = $2; }
  : Valve Switching Idle { $0 = SWS_VALVE_SWITCHING_IDLE; }
  : Overnight Flow Idle { $0 = SWS_OVERNIGHT_FLOW_IDLE; }
  : Overnight Flow Activate { $0 = SWS_OVERNIGHT_FLOW_ACTIVATE; }
  : Overnight Flow Execute { $0 = SWS_OVERNIGHT_FLOW_EXECUTE; }
  : Blank Steps Idle { $0 = SWS_BLANK_STEPS_IDLE; }
  : Blank Steps Activate { $0 = SWS_BLANK_STEPS_ACTIVATE; }
  : LC Steps Idle { $0 = SWS_LC_STEPS_IDLE; }
  : LC Steps Activate { $0 = SWS_LC_STEPS_ACTIVATE; }
  : Time Warp { $0 = SWS_TIMEWARP; }
  : Shutdown { $0 = SWS_SHUTDOWN; }
  ;
&Flag <unsigned char>
  : Set %d { $0 = $2; }
  ;
&Flag2 <unsigned char>
  : Set %d { $0 = $2; }
  ;
