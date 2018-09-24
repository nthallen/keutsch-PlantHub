#ifndef LICOR_H_INCLUDED
#define LICOR_H_INCLUDED

typedef struct __attribute__((__packed__)) {
  float CO2_mV;
  float CO2_ppm;
  float H2O_mV;
  float H2O_ppth;
  float Temp_mV;
  float P_kPa;
  uint16_t Status;
} licor_tm_t;

#define LICOR_FRESH 0x0001  // At least one record
#define LICOR_VFRESH 0x0002 // More than one record
#endif
