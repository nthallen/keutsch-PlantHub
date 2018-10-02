#include <string.h>
#include "LICOR_int.h"
#include "oui.h"
#include "nortlib.h"

const char *licor_path = "/dev/serusb1";
const char *licor_name = "LICOR";

int main( int argc, char **argv) {
  oui_init_options(argc, argv);
  nl_error(0, "Starting V1.0");
  { Selector S;
    Cmd_Selectee Cmd;
    licor_tm_t TMdata;
    LICOR_t LICOR(licor_path, &TMdata);
    TM_Selectee TM(licor_name, &TMdata, sizeof(licor_tm_t));
    S.add_child(&Cmd);
    S.add_child(&LICOR);
    S.add_child(&TM);
    S.event_loop();
  }
  nl_error(0, "Terminating");
}

LICOR_t::LICOR_t(const char *ser_dev, licor_tm_t *TMdata)
  : Ser_Sel(ser_dev, O_RDONLY|O_NONBLOCK, 58) {
  setup(4800, 8, 'n', 2, 10, 1);
  this->TMdata = TMdata;
  memset(TMdata, 0, sizeof(licor_tm_t));
  flush_input();
  flags |= TMgflag | Selector::Sel_Timeout;
  TO.Set(2, 0);
}

LICOR_t::~LICOR_t() {}

const int LICOR_t::TMgflag = Selector::gflag(0);

int LICOR_t::not_Lfloat(float &val) {
  for (; cp < nc && isspace(buf[cp]); ++cp) {
    if (buf[cp] == '\r' || buf[cp] == '\n') return 1;
  }
  return not_float(val);
}

int LICOR_t::ProcessData(int flag) {
  if (flag & TMgflag) {
    TMdata->Status &= ~(LICOR_FRESH | LICOR_VFRESH);
  }
  if (flag & (Selector::Sel_Read | Selector::Sel_Timeout)) {
    float CO2_mV;
    float CO2_ppm;
    float H2O_mV;
    float H2O_ppth;
    float Temp_mV;
    float P_kPa;
    if (fillbuf()) return 1;
    cp = 0;
    if (not_Lfloat(CO2_mV) ||
        not_Lfloat(CO2_ppm) ||
        not_Lfloat(H2O_mV) ||
        not_Lfloat(H2O_ppth) ||
        not_Lfloat(Temp_mV) ||
        not_Lfloat(P_kPa) ||
        not_str("\r\n")) {
      if (cp >= nc) {
        if (TO.Expired()) {
          report_err("Timeout reading from LICOR");
        } else {
          return 0;
        }
      } else {
        cp = 0;
        if (not_found('\r')) return 0;
        if (cp < nc && buf[cp] == '\n') ++cp;
      }
      consume(cp);
      TO.Set(2,0);
    } else {
      TMdata->CO2_mV = CO2_mV;
      TMdata->CO2_ppm = CO2_ppm;
      TMdata->H2O_mV = H2O_mV;
      TMdata->H2O_ppth = H2O_ppth;
      TMdata->Temp_mV = Temp_mV;
      TMdata->P_kPa = P_kPa;
      if (TMdata->Status | LICOR_FRESH) {
        TMdata->Status |= LICOR_VFRESH;
      }
      TMdata->Status |= LICOR_FRESH;
      report_ok();
      consume(cp);
      TO.Set(2,0);
    }
  }
  return 0;
}

Timeout *LICOR_t::GetTimeout() {
  return &TO;
}

