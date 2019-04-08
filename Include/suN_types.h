#ifdef GAUGE_SUN
  #ifdef REPR_ADJOINT
    #if NG == 2
      #include "su2adj_types.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 2
      #include "su2_types.h"
    #elif NG == 3
      #include "su3_types.h"
    #endif
  #elif defined REPR_SYMMETRIC
    #if NG == 3
      #include "su3as_types.h"
    #endif
  #endif
#elif defined GAUGE_SPN
  #ifdef REPR_ADJOINT
    #if NG == 4
      #include "sp4adj_types.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 4
      #include "sp4_types.h"
    #endif
  #endif
#endif