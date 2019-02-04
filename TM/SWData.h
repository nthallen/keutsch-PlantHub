/* SWData.h */
#ifndef SWDATA_H_INCLUDED
#define SWDATA_H_INCLUDED

typedef struct __attribute__((__packed__)) {
  unsigned char SWStat;
  unsigned char Flag;
} SWData_t;
extern SWData_t SWData;

#define SWS_VALVE_SWITCHING_ACTIVE 1
#define SWS_VALVE_SWITCHING_IDLE 2
#define SWS_AUTO_START_FLOW_IDLE 3
#define SWS_AUTO_START_FLOW_ACTIVE 4
#define SWS_AUTO_START_FLOW_EXECUTE 5
#define SWS_FCC_TEST_IDLE 6
#define SWS_FCC_TEST 7
#define SWS_CO2_CONTROL_IDLE 8
#define SWS_CO2_CONTROL_INITIALIZE 9
#define SWS_TIMEWARP 254
#define SWS_SHUTDOWN 255

#endif
