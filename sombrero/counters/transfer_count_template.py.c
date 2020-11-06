#include "fm_defs.h"
#include "memory_base.py.h"

#ifndef MKPYMOD
#include "transfer_count.h"
#ifndef TESTCOUNTERS
#include "libhr_defines_interface.h" // provides T,X,Y,Z and NF
#endif
#endif

_FD(spinor_field_memory_transfer, halo_sites() * spinor_size());

_FD(Dphi_memory_transfer, spinor_field_memory_transfer());

#ifdef MKPYMOD
def g5Cphi_eopre_sq_memory_transfer():
#else
float g5Cphi_eopre_sq_memory_transfer(){
#endif
    return 2 * 2 * Dphi_memory_transfer();
#ifndef MKPYMOD
}
#endif

// This complicated pattern is followed for uniformity with the
// flop and the memory counters.
// Neglecting communication of scalar values.
#ifdef MKPYMOD
def cg_out_of_loop_memory_transfer(operator_memory_transfer):
#else
float cg_out_of_loop_memory_transfer(float operator_memory_transfer){
#endif
    return operator_memory_transfer;
#ifndef MKPYMOD
}
#endif

// This complicated pattern is followed for uniformity with the
// flop and the memory counters.
// Neglecting communication of scalar values.
#ifdef MKPYMOD
def cg_iteration_memory_transfer(operator_memory_transfer):
#else
float cg_iteration_memory_transfer(float operator_memory_transfer) {
#endif
    return operator_memory_transfer;
#ifndef MKPYMOD
}
#endif
