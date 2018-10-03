#ifndef LICOR_H_INCLUDED
#define LICOR_H_INCLUDED
#include <stdint.h>

typedef struct __attribute__((__packed__)) {
  float CO2_mV;
  float CO2_ppm;
  float H2O_mV;
  float H2O_ppth;
  float Temp_mV;
  float P_kPa;
  uint16_t Status;
} licor_tm_t;

// LICOR_FRESH means At least one record
// LICOR_VFRESH means More than one record
#define LICOR_FRESH 1
#define LICOR_VFRESH 2
#endif
