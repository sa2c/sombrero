#ifndef __MEMORY_COUNT_H_
#define __MEMORY_COUNT_H_

float local_apply_BCs_on_spinor_field_memory();
float local_g5Cphi_eopre_sq_memory(float bc_operation_memory);
float local_cg_iteration_memory(float local_operator_memory);
float local_cg_out_of_loop_memory(float local_operator_memory);
float cg_local_GB_used(int real_iterations);

#endif // __MEMORY_COUNT_H_
