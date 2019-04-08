#ifdef GAUGE_SUN
  #ifdef REPR_ADJOINT
    #if NG == 2
      #include "su2adj_repr_func.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 2
      #include "su2_repr_func.h"
    #elif NG == 3
      #include "su3_repr_func.h"
    #endif
  #elif defined REPR_SYMMETRIC
    #if NG == 3
      #include "su3as_repr_func.h"
    #endif
  #endif
#elif defined GAUGE_SPN
  #ifdef REPR_ADJOINT
    #if NG == 4
      #include "sp4adj_repr_func.h"
    #endif
  #elif defined REPR_FUNDAMENTAL
    #if NG == 4
      #include "sp4_repr_func.h"
    #endif
  #endif
#endif