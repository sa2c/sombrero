#ifdef GAUGE_SUN
  #ifdef REPR_ADJOINT
    #if NG == 2
      #include "su2adj_exp.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 2
      #include "su2_exp.h"
    #elif NG == 3
      #include "su3_exp.h"
    #endif
  #elif defined REPR_SYMMETRIC
    #if NG == 3
      #include "su3as_exp.h"
    #endif
  #endif
#elif defined GAUGE_SPN
  #ifdef REPR_ADJOINT
    #if NG == 4
      #include "sp4adj_exp.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 4
      #include "sp4_exp.h"
    #endif
  #endif
#endif