#include "nl.h"
#include "SB.h"

subbuspp *FCC0;

void FCC0_init() {
  FCC0 = new subbuspp("SB", "serusb");
  if (FCC0->load()) {
    uint16_t BdID = FCC0->sbrd(3);
    if (BdID != PLANT_FCC_BDID) {
      msg(2, "Incorrect BdID %d: expected %d", BdID, PLANT_FCC_BDID);
    }
  } else {
    msg(2, "FCC Subbus library not found");
    delete FCC0;
    FCC0 = 0;
  }
}
