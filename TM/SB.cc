#include "nortlib.h"
#include "SB.h"

subbuspp *FCC0;

void FCC0_init() {
  FCC0 = new subbuspp("SB");
  if (FCC0->load()) {
    uint16_t BdID = FCC0->sbrd(2);
    uint16_t InstID = FCC0->sbrd(5);
    if (BdID != PLANT_FCC_BDID || InstID != PLANT_INST_ID) {
      nl_error(2, "Incorrect BdID:InstID %d:%d expected %d:%d",
        BdID, InstID, PLANT_FCC_BDID, PLANT_INST_ID);
    }
  } else {
    nl_error(2, "FCC Subbus library not found");
    delete FCC0;
    FCC0 = 0;
  }
}
