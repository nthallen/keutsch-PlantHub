%{
  #include <stdint.h>
  #include "SB.h"
  /* fcc.cmd */
  #ifdef SERVER
    int SCCM_Span[3] = { 10, 500, 5000};
  #endif
%}

&command
  : Flow &flowchan Set Volts %f (Enter value in Volts) * {
      double dsp = $5*65536./5.;
      uint16_t usp;
      
      if (dsp > 65535) usp = 65535;
      else if (dsp < 0) usp = 0;
      else usp = dsp;
      if (FCC0) FCC0->sbwr(0x14+$2, usp);
      else nl_error(2, "FCC not present");
    }
  : Flow &flowchan Set SCCM %f (Enter value in SCCM) * {
      double dsp;
      uint16_t usp;
      
      if ($2 >= 0 && $2 <= 2) {
        dsp = $5*65536./SCCM_Span[$2];
        if (dsp > 65535) usp = 65535;
        else if (dsp < 0) usp = 0;
        else usp = dsp;
        if (FCC0) FCC0->sbwr(0x14+$2, usp);
        else nl_error(2, "FCC not present");
      } else nl_error(2, "Invalid DAC Channel: %d", $2);
    }
  : Flow &flowchan Solenoid &OpenCloseCtrl * {
      uint16_t cmd = ($2<<2) + $4;
      nl_error(0, "FCC command is %u", cmd);
      if (FCC0) FCC0->sbwr(0x18, cmd);
      else nl_error(2, "FCC not present");
    }
  : Valve &ValveSelect &OpenClose * {
      uint16_t cmd = 38 + $2*2 + $3;
      if (FCC0) FCC0->sbwr(0x18, cmd);
      else nl_error(2, "FCC not present");
    }
  ;

&flowchan <uint16_t>
  : CO2 { $0 = 0; }
  : OVOC { $0 = 1; }
  : Zero { $0 = 2; }
  ;

&OpenCloseCtrl <uint16_t>
  : Open { $0 = 2; }
  : Close { $0 = 1; }
  : Control { $0 = 0; }
  ;

&ValveSelect <uint16_t>
  : 1 { $0 = 0; }
  : 2 { $0 = 1; }
  ;

&OpenClose <uint16_t>
  : On { $0 = 1; }
  : Off { $0 = 0; }
  ;
