/* SWData.h */
#ifndef SWDATA_H_INCLUDED
#define SWDATA_H_INCLUDED

typedef struct __attribute__((__packed__)) {
  unsigned char SWStat;
  unsigned char Flag;
} SWData_t;
extern SWData_t SWData;

#define SWS_PCTRL_IDLE 1
#define SWS_PCTRL_ACTIVE 2
#define SWS_VALVE_SWITCHING_ACTIVE 3
#define SWS_VALVE_SWITCHING_IDLE 4
#define SWS_TIMEWARP 254
#define SWS_SHUTDOWN 255

#endif
