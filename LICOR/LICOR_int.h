#include LICOR_INT_H_INCLUDED
#define LICOR_INT_H_INCLUDED

extern const char *licor_path;
extern const char *licor_name;

#ifdef __cplusplus

#include "SerSelector.h"
#include "LICOR.h"
#include "collect.h"

class LICOR_t : public Ser_Sel {
  public:
    LICOR_t(const char *ser_dev, licor_tm_t *data);
    ~LICOR_t();
    int ProcessData(int flag);
    Timeout *GetTimeout();
    const int TMgflag = Selector::gflag(0);
  private:
    Timeout TO;
    licor_tm_t *TMdata;
};

#endif /* __cplusplus */

#endif /* LICOR_INT_H_INCLUDED */
