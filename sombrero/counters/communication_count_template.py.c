#include "fm_defs.h"
#include "memory_base.py.h"

#ifndef MKPYMOD
#include "communication_count.h"
#include "libhr_defines_interface.h" // provides T,X,Y,Z, NF and border sizes
#endif

_FD(spinor_field_communication, 2 * halo_sites_e() * spinor_size());

_FD(Dphi_communication, spinor_field_communication());

#ifdef MKPYMOD
def g5Cphi_eopre_sq_communication():
#else
float g5Cphi_eopre_sq_communication(){
#endif
    return 2 * 2 * Dphi_communication();
#ifndef MKPYMOD
}
#endif

// This complicated pattern is followed for uniformity with the
// flop and the memory counters.
// Neglecting communication of scalar values.
#ifdef MKPYMOD
def cg_out_of_loop_communication(operator_communication):
#else
float cg_out_of_loop_communication(float operator_communication){
#endif
    return operator_communication;
#ifndef MKPYMOD
}
#endif

// This complicated pattern is followed for uniformity with the
// flop and the memory counters.
// Neglecting communication of scalar values.
#ifdef MKPYMOD
def cg_iteration_communication(operator_communication):
#else
float cg_iteration_communication(float operator_communication) {
#endif
    return operator_communication;
#ifndef MKPYMOD
}
#endif
