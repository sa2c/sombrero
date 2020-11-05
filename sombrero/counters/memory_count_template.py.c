#include "fc_defs.h"
#include "memory_base.py.h"

#ifndef MKPYMOD
#include "memory_count.h"
#ifndef TESTING
#include "global.h"
#endif
#endif

_FD(ldl_t_size, 2 * NF * (2 * NF + 1) * complex_size());
_FD(gauge_matrix_size, NF *vector_size());

_FD(local_spinor_field_e_memory,
    spinor_size() * (local_sites_e() + halo_sites_e()));

_FD(local_spinor_field_e_nohalo_memory, spinor_size() * local_sites_e());

// the only gauge links needed from the halo are 1 per halo site
// only for the halos in the "back" direction.
_FD(local_gauge_field_memory,
    gauge_matrix_size() * (local_sites() * 4 + halo_sites() / 2));

_FD(local_clover_field_e_memory, gauge_matrix_size() * local_sites_e() * 4);

_FD(local_ldl_field_e_memory, ldl_t_size() * local_sites_e());

// void Dphi_(spinor_field *out, spinor_field *in)
// gauge field implicitly
_FD(local_Dphi_memory,
    2 * local_spinor_field_e_memory() + local_gauge_field_memory());

// static void Cphi_(double mass, spinor_field *dptr, spinor_field *sptr, int
// assign) cl_term implicit
_FD(local_Cphi_memory,
    2 * local_spinor_field_e_memory() + local_clover_field_e_memory());

// static void Cphi_inv_(double mass, spinor_field *dptr, spinor_field *sptr,
// int assign) cl_ldl field implicitly
_FD(local_Cphi_inv_memory,
    2 * local_spinor_field_e_memory() + local_ldl_field_e_memory());

// TO UNDERSTAND - boundary conditions
_FD(tslice_ext_sites, (X + 2) * (Y + 2) * (Z + 2));

_FD(local_spinor_field_e_tslice_memory, spinor_size() * tslice_ext_sites());

#ifdef MKPYMOD
def local_apply_BCs_on_spinor_field_memory():
#else
float local_apply_BCs_on_spinor_field_memory(){
#endif
   return local_spinor_field_e_tslice_memory();
#ifndef MKPYMOD
}
#endif

_FD(local_spinor_field_minus_f_memory, 2 * local_spinor_field_e_memory())

// void Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr)
// Assuming there is no data reuse between functions
#ifdef MKPYMOD
def local_Cphi_eopre_memory(bc_operation_memory):
#else
float local_Cphi_eopre_memory(float bc_operation_memory) {
#endif
      return (bc_operation_memory + local_Dphi_memory() +
              local_Cphi_inv_memory() + bc_operation_memory +
              local_Dphi_memory() + local_spinor_field_minus_f_memory() +
              local_Cphi_memory() + bc_operation_memory);
#ifndef MKPYMOD
}
#endif

_FD(local_spinor_field_g5_assign_f_memory, 2 * local_spinor_field_e_memory())

// Assuming there is no data reuse between functions
#ifdef MKPYMOD
def local_g5Cphi_eopre_sq_memory(bc_operation_memory):
#else
float local_g5Cphi_eopre_sq_memory(float bc_operation_memory) {
#endif
    return (local_Cphi_eopre_memory(bc_operation_memory) +
            local_spinor_field_g5_assign_f_memory());
#ifndef MKPYMOD
}
#endif
// Assuming there is no data reuse between functions
//
//
_FD(local_spinor_field_sqnorm_f_memory,
    REAL_SIZE + local_spinor_field_e_nohalo_memory());

_FD(local_spinor_field_mul_add_assign_f_memory,
    REAL_SIZE + 2 * local_spinor_field_e_memory());

_FD(local_spinor_field_sub_f_memory, 2 * local_spinor_field_e_memory());

_FD(local_spinor_field_copy_f_memory, 2 * local_spinor_field_e_memory());

#ifdef MKPYMOD
def local_cg_out_of_loop_memory(local_operator_memory):
#else
float local_cg_out_of_loop_memory(float local_operator_memory) {
#endif
    return (local_spinor_field_sqnorm_f_memory() + local_operator_memory +
              local_spinor_field_mul_add_assign_f_memory() +
              local_spinor_field_sub_f_memory() +
              local_spinor_field_copy_f_memory() +
              local_spinor_field_mul_add_assign_f_memory() +
              local_spinor_field_sub_f_memory() +
              local_spinor_field_sqnorm_f_memory());
#ifndef MKPYMOD
}
#endif

_FD(local_spinor_prod_re_f_memory,
    REAL_SIZE + 2 * local_spinor_field_e_memory());
_FD(local_spinor_field_mul_f_memory,
    REAL_SIZE + 2 * local_spinor_field_e_memory());

#ifdef MKPYMOD
def local_cg_iteration_memory(local_operator_memory):
#else
float local_cg_iteration_memory(float local_operator_memory) {
#endif
    return (local_operator_memory + local_spinor_prod_re_f_memory() +
            local_spinor_field_mul_add_assign_f_memory() +
            local_spinor_field_mul_add_assign_f_memory() +
            local_spinor_field_sqnorm_f_memory() +
            local_spinor_field_mul_f_memory() +
            local_spinor_field_mul_add_assign_f_memory());
#ifndef MKPYMOD
}
#endif
