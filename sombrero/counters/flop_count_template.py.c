#include "fc_defs.h"

#ifndef PYTHON
#include "flop_count.h"
#endif

#ifdef PYTHON
NF = 4 // just to have it declared.
#endif

#ifdef REPR_ADJOINT
    _FD(matrix_mul_flops, 4 * NF * NF - 2 * NF);
#else
_FD(matrix_mul_flops, 8 * NF * NF - 2 * NF);
#endif

_FD(real_op_flop, 1);
_FD(complex_sum_flops, 2);
_FD(complex_mul_flops, 6);
_FD(complex_mulR_flops, 2); // capital R for legibility (hopefully);
_FD(complex_norm_flops, 3); // only re ;

// single vector
_FD(vector_sum_flops, complex_sum_flops() * NF);
_FD(vector_mul_flops, complex_mulR_flops() * NF); // mul uses mulR ;
_FD(vector_norm_flops, complex_norm_flops() * NF +
                       real_op_flop() * NF); // + NF sums (not NF-1)

// single spinor
_FD(spinor_norm_flops, 4 * vector_norm_flops());
_FD(spinor_sum_flops, 4 * vector_sum_flops());
_FD(spinor_mul_flops, 4 * vector_mul_flops()); // uses mulR
_FD(spinor_mul_add_assign_flops, spinor_mul_flops() +
                                 spinor_sum_flops());

// spinor field related flop count, for a single site
_FD(site_spinor_field_sqnorm_f_flops,
    1.0 / 2 * spinor_norm_flops()); // mult + add
_FD(site_spinor_field_mul_add_assign_f_flops,
    1.0 / 2 * spinor_mul_add_assign_flops());
_FD(site_spinor_field_mul_f_flops, 1.0 / 2 * spinor_mul_flops());
_FD(site_spinor_field_sub_f_flops, 1.0 / 2 * spinor_sum_flops());
_FD(site_spinor_field_prod_re_flop, site_spinor_field_sqnorm_f_flops());
_FD(site_spinor_field_add_assign_f_flops, 1.0 / 2 * spinor_sum_flops());
_FD(site_spinor_field_minus_f_flops, 1.0 / 2 * spinor_mul_flops()); // as a mul

// see Dphi.c
_FD(site_Dphi_flops, 1.0 / 2 * (16 * matrix_mul_flops() +
                                45 * vector_sum_flops()));

_FD(site_Cphi_flops, 1.0 / 2 * (8 * matrix_mul_flops() +
                                4 * vector_sum_flops() +
                                spinor_mul_add_assign_flops()));

_FD(site_Cphi_assign_flops, 1.0 / 2 * (8 * matrix_mul_flops() +
                                       4 * vector_sum_flops() +
                                       spinor_mul_add_assign_flops() +
                                       spinor_sum_flops()));

_FD(clover_matsize, 2 * NF);

_FD(forward_substitution_loop_flop,
    clover_matsize() * (clover_matsize() - 1) / 2 *
    ((3 * real_op_flop()) + // n = i*(i+1)/2+i)
     (complex_mul_flops() + complex_sum_flops()) + // _complex_mul_sub_assign();
     (complex_mul_flops() + complex_sum_flops()))); // _complex_mul_sub_assign();

#ifdef PYTHON
def backward_substitution_loop_flop():
    res = 0
#else
float backward_substitution_loop_flop() {
    int res = 0;
#endif

#ifdef PYTHON
    for i in range(clover_matsize()-1,-1,-1):
#else
    for (int i = clover_matsize() - 1; i >= 0; i--) {
#endif
        res += 3; // n = i*(i+1)/2+i);
        res += complex_mulR_flops(); //_complex_mulr());
        res += complex_mulR_flops(); //_complex_mulr());
#ifdef PYTHON
        for k in range(i + 1, N()):
#else
        for (int k = i + 1; k < clover_matsize(); k++) {
#endif
            res += 3; // n = k*(k+1)/2+i);
            // _complex_mul_sub_assign());
            res += complex_mul_flops() + complex_sum_flops();
            // _complex_mul_sub_assign());
            res += complex_mul_flops() + complex_sum_flops();

#ifndef PYTHON
        }
    }
#endif
    return res;
#ifndef PYTHON
}
#endif

_FD(site_Cphi_inv_flops, 1.0 / 2 *
    (forward_substitution_loop_flop() +
     backward_substitution_loop_flop()));

_FD(site_Cphi_inv_assign_flops,
    1.0 / 2 *
    (forward_substitution_loop_flop() +
     backward_substitution_loop_flop() + spinor_sum_flops()));

#ifdef PYTHON
def site_g5Cphi_eopre_sq_flops():
#else
float site_g5Cphi_eopre_sq_flops() {
    float site_g5Cphi_eopre_flops;
#endif
    site_g5Cphi_eopre_flops = (site_Dphi_flops() +
                               site_Cphi_inv_flops() +
                               site_Dphi_flops() +
                               site_spinor_field_minus_f_flops() +
                               site_Cphi_assign_flops());
    return 2 * site_g5Cphi_eopre_flops;
#ifndef PYTHON
}
#endif

#ifdef PYTHON
def cg_out_of_loop_flops_per_site(site_operator_flops):
#else
float cg_out_of_loop_flops_per_site(float site_operator_flops) {
#endif
    return (site_spinor_field_sqnorm_f_flops() + site_operator_flops +
            site_spinor_field_mul_add_assign_f_flops() +
            site_spinor_field_sub_f_flops() +
            site_spinor_field_sqnorm_f_flops() + site_operator_flops +
            site_spinor_field_mul_add_assign_f_flops() +
            site_spinor_field_sub_f_flops() +
            2 * site_spinor_field_sqnorm_f_flops());
#ifndef PYTHON
}
#endif

#ifdef PYTHON
def cg_iteration_flops_per_site(site_operator_flops):
#else
float cg_iteration_flops_per_site(float site_operator_flops) {
#endif
    return (site_operator_flops + site_spinor_field_prod_re_flop() +
            site_spinor_field_mul_add_assign_f_flops() +
            site_spinor_field_mul_add_assign_f_flops() +
            site_spinor_field_sqnorm_f_flops() +
            site_spinor_field_mul_f_flops() +
            site_spinor_field_add_assign_f_flops() +
            site_spinor_field_mul_f_flops() +
            site_spinor_field_add_assign_f_flops());
#ifndef PYTHON
}
#endif
