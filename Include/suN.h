#ifdef GAUGE_SUN
  #ifdef REPR_ADJOINT
    #if NG == 2
      #include "su2adj.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 2
      #include "su2.h"
    #elif NG == 3
      #include "su3.h"
    #endif
  #elif defined REPR_SYMMETRIC
    #if NG == 3
      #include "su3as.h"
    #endif
  #endif
#elif defined GAUGE_SPN
  #ifdef REPR_ADJOINT
    #if NG == 4
      #include "sp4adj.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 4
      #include "sp4.h"
    #endif
  #endif
#endif