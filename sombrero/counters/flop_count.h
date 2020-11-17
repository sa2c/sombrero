#ifndef __FLOP_COUNT_H_
#define __FLOP_COUNT_H_

float site_g5Cphi_eopre_sq_flops();
float cg_iteration_flops_per_site(float site_operator_flops);
float cg_out_of_loop_flops_per_site(float site_operator_flops);
float cg_Gflops(float real_iterations);

#endif // __FLOP_COUNT_H_
