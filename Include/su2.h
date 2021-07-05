#define _vector_zero_g(r)\
 _complex_0((r).c[0]);\
 _complex_0((r).c[1])

#define _suNg_zero(u)\
 _complex_0((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3])

#define _algebra_vector_zero_g(r)\
 (r).c[0]=0.;\
 (r).c[1]=0.;\
 (r).c[2]=0.

#define _suNg_unit(u)\
 _complex_1((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_1((u).c[3])

#define _vector_zero_f(r)\
 _complex_0((r).c[0]);\
 _complex_0((r).c[1])

#define _vector_minus_f(r,s)\
 _complex_minus((r).c[0],(s).c[0]);\
 _complex_minus((r).c[1],(s).c[1])

#define _vector_i_plus_f(r,s)\
 _complex_i_plus((r).c[0],(s).c[0]);\
 _complex_i_plus((r).c[1],(s).c[1])

#define _vector_i_minus_f(r,s)\
 _complex_i_minus((r).c[0],(s).c[0]);\
 _complex_i_minus((r).c[1],(s).c[1])

#define _vector_mul_f(r,k,s)\
 _complex_mulr((r).c[0],(k),(s).c[0]);\
 _complex_mulr((r).c[1],(k),(s).c[1])

#define _vector_mulc_f(r,z,s)\
 _complex_mul((r).c[0],(z),(s).c[0]);\
 _complex_mul((r).c[1],(z),(s).c[1])

#define _vector_mulc_star_f(r,z,s)\
 _complex_mul_star((r).c[0],(s).c[0],(z));\
 _complex_mul_star((r).c[1],(s).c[1],(z))

#define _vector_add_f(r,s1,s2)\
 _complex_add((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_add((r).c[1],(s1).c[1],(s2).c[1])

#define _vector_sub_f(r,s1,s2)\
 _complex_sub((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_sub((r).c[1],(s1).c[1],(s2).c[1])

#define _vector_i_add_f(r,s1,s2)\
 _complex_i_add((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_i_add((r).c[1],(s1).c[1],(s2).c[1])

#define _vector_i_sub_f(r,s1,s2)\
 _complex_i_sub((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_i_sub((r).c[1],(s1).c[1],(s2).c[1])

#define _vector_add_assign_f(r,s)\
 _complex_add_assign((r).c[0],(s).c[0]);\
 _complex_add_assign((r).c[1],(s).c[1])

#define _vector_sub_assign_f(r,s)\
 _complex_sub_assign((r).c[0],(s).c[0]);\
 _complex_sub_assign((r).c[1],(s).c[1])

#define _vector_i_add_assign_f(r,s)\
 _complex_i_add_assign((r).c[0],(s).c[0]);\
 _complex_i_add_assign((r).c[1],(s).c[1])

#define _vector_i_sub_assign_f(r,s)\
 _complex_i_sub_assign((r).c[0],(s).c[0]);\
 _complex_i_sub_assign((r).c[1],(s).c[1])

#define _vector_prod_re_f(k,r,s)\
 (k)=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_re((r).c[1],(s).c[1])

#define _vector_prod_im_f(k,r,s)\
 (k)=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_im((r).c[1],(s).c[1])

#define _vector_mulc_add_assign_f(r,z,s)\
 _complex_mul_assign((r).c[0],(z),(s).c[0]);\
 _complex_mul_assign((r).c[1],(z),(s).c[1])

#define _vector_mul_add_assign_f(r,k,s)\
 _complex_mulr_assign((r).c[0],(k),(s).c[0]);\
 _complex_mulr_assign((r).c[1],(k),(s).c[1])

#define _vector_lc_f(r,k1,s1,k2,s2)\
 _complex_rlc((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _complex_rlc((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1])

#define _vector_lc_add_assign_f(r,k1,s1,k2,s2)\
 _complex_rlc_assign((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _complex_rlc_assign((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1])

#define _vector_clc_f(r,z1,s1,z2,s2)\
 _complex_clc((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _complex_clc((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1])

#define _vector_clc_add_assign_f(r,z1,s1,z2,s2)\
 _complex_clc_assign((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _complex_clc_assign((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1])

#define _vector_prod_assign_f(z,r,s)\
 _complex_prod_assign((z),(r).c[0],(s).c[0]);\
 _complex_prod_assign((z),(r).c[1],(s).c[1])

#define _vector_prod_add_assign_re_f(k,r,s)\
 (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_re((r).c[1],(s).c[1])

#define _vector_prod_add_assign_im_f(k,r,s)\
 (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_im((r).c[1],(s).c[1])

#define _vector_prod_sub_assign_re_f(k,r,s)\
 (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)-=_complex_prod_re((r).c[1],(s).c[1])

#define _vector_prod_sub_assign_im_f(k,r,s)\
 (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)-=_complex_prod_im((r).c[1],(s).c[1])

#define _suNf_multiply(r,u,s)\
 _complex_mul((r).c[0],(u).c[0],(s).c[0]);\
 _complex_mul_assign((r).c[0],(u).c[1],(s).c[1]);\
 _complex_mul((r).c[1],(u).c[2],(s).c[0]);\
 _complex_mul_assign((r).c[1],(u).c[3],(s).c[1])

#define _suNf_inverse_multiply(r,u,s)\
 _complex_mul_star((r).c[0],(s).c[0],(u).c[0]);\
 _complex_mul_star_assign((r).c[0],(s).c[1],(u).c[2]);\
 _complex_mul_star((r).c[1],(s).c[0],(u).c[1]);\
 _complex_mul_star_assign((r).c[1],(s).c[1],(u).c[3])

#define _suNf_zero(u)\
 _complex_0((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3])

#define _suNf_dagger(u,v)\
 _complex_star((u).c[0],(v).c[0]);\
 _complex_star((u).c[1],(v).c[2]);\
 _complex_star((u).c[2],(v).c[1]);\
 _complex_star((u).c[3],(v).c[3])

#define _suNf_times_suNf(u,v,w)\
 _complex_mul((u).c[0],(v).c[0],(w).c[0]);\
 _complex_mul_assign((u).c[0],(v).c[1],(w).c[2]);\
 _complex_mul((u).c[1],(v).c[0],(w).c[1]);\
 _complex_mul_assign((u).c[1],(v).c[1],(w).c[3]);\
 _complex_mul((u).c[2],(v).c[2],(w).c[0]);\
 _complex_mul_assign((u).c[2],(v).c[3],(w).c[2]);\
 _complex_mul((u).c[3],(v).c[2],(w).c[1]);\
 _complex_mul_assign((u).c[3],(v).c[3],(w).c[3])

#define _suNf_times_suNf_dagger(u,v,w)\
 _complex_mul_star((u).c[0],(v).c[0],(w).c[0]);\
 _complex_mul_star_assign((u).c[0],(v).c[1],(w).c[1]);\
 _complex_mul_star((u).c[1],(v).c[0],(w).c[2]);\
 _complex_mul_star_assign((u).c[1],(v).c[1],(w).c[3]);\
 _complex_mul_star((u).c[2],(v).c[2],(w).c[0]);\
 _complex_mul_star_assign((u).c[2],(v).c[3],(w).c[1]);\
 _complex_mul_star((u).c[3],(v).c[2],(w).c[2]);\
 _complex_mul_star_assign((u).c[3],(v).c[3],(w).c[3])

#define _suNf_dagger_times_suNf(u,v,w)\
 _complex_mul_star((u).c[0],(w).c[0],(v).c[0]);\
 _complex_mul_star_assign((u).c[0],(w).c[2],(v).c[2]);\
 _complex_mul_star((u).c[1],(w).c[1],(v).c[0]);\
 _complex_mul_star_assign((u).c[1],(w).c[3],(v).c[2]);\
 _complex_mul_star((u).c[2],(w).c[0],(v).c[1]);\
 _complex_mul_star_assign((u).c[2],(w).c[2],(v).c[3]);\
 _complex_mul_star((u).c[3],(w).c[1],(v).c[1]);\
 _complex_mul_star_assign((u).c[3],(w).c[3],(v).c[3])

#define _suNf_minus(u,v)\
 _complex_minus((u).c[0],(v).c[0]);\
 _complex_minus((u).c[1],(v).c[1]);\
 _complex_minus((u).c[2],(v).c[2]);\
 _complex_minus((u).c[3],(v).c[3])

#define _suNf_mul(u,r,v)\
 _complex_mulr((u).c[0],(r),(v).c[0]);\
 _complex_mulr((u).c[1],(r),(v).c[1]);\
 _complex_mulr((u).c[2],(r),(v).c[2]);\
 _complex_mulr((u).c[3],(r),(v).c[3])

#define _suNf_mul_assign(u,r)\
 (u).c[0]*=(r);\
 (u).c[1]*=(r);\
 (u).c[2]*=(r);\
 (u).c[3]*=(r)

#define _suNf_mulc(u,r,v)\
 _complex_mul((u).c[0],(r),(v).c[0]);\
 _complex_mul((u).c[1],(r),(v).c[1]);\
 _complex_mul((u).c[2],(r),(v).c[2]);\
 _complex_mul((u).c[3],(r),(v).c[3])

#define _suNf_add_assign(u,v)\
 _complex_add_assign((u).c[0],(v).c[0]);\
 _complex_add_assign((u).c[1],(v).c[1]);\
 _complex_add_assign((u).c[2],(v).c[2]);\
 _complex_add_assign((u).c[3],(v).c[3])

#define _suNfc_multiply(a,b,c)\
 _suNf_multiply(a,b,c)

#define _suNfc_inverse_multiply(a,b,c)\
 _suNf_inverse_multiply(a,b,c)

#define _suNfc_zero(a)\
 _suNf_zero(a)

#define _suNfc_dagger(u,v)\
 _suNf_dagger(u,v)

#define _suNfc_mul(u,r,v)\
 _suNf_mul(u,r,v)

#define _suNfc_mul_assign(u,r)\
 _suNf_mul_assign(u,r)

#define _suNffull_multiply(a,b,c)\
 _suNf_multiply(a,b,c)

#define _suNffull_mul(a,b,c)\
 _suNf_mul(a,b,c)

#define _suNffull_inverse_multiply(a,b,c)\
 _suNf_inverse_multiply(a,b,c)

#define _suNf_expand(a,b)\
 _suNf_copy(a,b)

#define _spinor_zero_f(r)\
 _vector_zero_f((r).c[0]);\
 _vector_zero_f((r).c[1]);\
 _vector_zero_f((r).c[2]);\
 _vector_zero_f((r).c[3])

#define _spinor_g5_f(s,r)\
 (s).c[0]=(r).c[0];\
 (s).c[1]=(r).c[1];\
 _vector_minus_f((s).c[2],(r).c[2]);\
 _vector_minus_f((s).c[3],(r).c[3])

#define _spinor_g5_assign_f(r)\
 _vector_minus_f((r).c[2],(r).c[2]);\
 _vector_minus_f((r).c[3],(r).c[3])

#define _spinor_minus_f(s,r)\
 _vector_minus_f((s).c[0],(r).c[0]);\
 _vector_minus_f((s).c[1],(r).c[1]);\
 _vector_minus_f((s).c[2],(r).c[2]);\
 _vector_minus_f((s).c[3],(r).c[3])

#define _spinor_mul_f(r,k,s)\
 _vector_mul_f((r).c[0],k,(s).c[0]);\
 _vector_mul_f((r).c[1],k,(s).c[1]);\
 _vector_mul_f((r).c[2],k,(s).c[2]);\
 _vector_mul_f((r).c[3],k,(s).c[3])

#define _spinor_mulc_f(r,z,s)\
 _vector_mulc_f((r).c[0],z,(s).c[0]);\
 _vector_mulc_f((r).c[1],z,(s).c[1]);\
 _vector_mulc_f((r).c[2],z,(s).c[2]);\
 _vector_mulc_f((r).c[3],z,(s).c[3])

#define _spinor_mulc_add_assign_f(r,z,s)\
 _vector_mulc_add_assign_f((r).c[0],(z),(s).c[0]);\
 _vector_mulc_add_assign_f((r).c[1],(z),(s).c[1]);\
 _vector_mulc_add_assign_f((r).c[2],(z),(s).c[2]);\
 _vector_mulc_add_assign_f((r).c[3],(z),(s).c[3])

#define _spinor_mul_add_assign_f(r,k,s)\
 _vector_mul_add_assign_f((r).c[0],(k),(s).c[0]);\
 _vector_mul_add_assign_f((r).c[1],(k),(s).c[1]);\
 _vector_mul_add_assign_f((r).c[2],(k),(s).c[2]);\
 _vector_mul_add_assign_f((r).c[3],(k),(s).c[3])

#define _spinor_lc_f(r,k1,s1,k2,s2)\
 _vector_lc_f((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _vector_lc_f((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]);\
 _vector_lc_f((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]);\
 _vector_lc_f((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

#define _spinor_lc_add_assign_f(r,k1,s1,k2,s2)\
 _vector_lc_add_assign_f((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _vector_lc_add_assign_f((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]);\
 _vector_lc_add_assign_f((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]);\
 _vector_lc_add_assign_f((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

#define _spinor_clc_f(r,z1,s1,z2,s2)\
 _vector_clc_f((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _vector_clc_f((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]);\
 _vector_clc_f((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]);\
 _vector_clc_f((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

#define _spinor_clc_add_assign_f(r,z1,s1,z2,s2)\
 _vector_clc_add_assign_f((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _vector_clc_add_assign_f((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]);\
 _vector_clc_add_assign_f((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]);\
 _vector_clc_add_assign_f((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

#define _spinor_add_f(r,s1,s2)\
 _vector_add_f((r).c[0],(s1).c[0],(s2).c[0]);\
 _vector_add_f((r).c[1],(s1).c[1],(s2).c[1]);\
 _vector_add_f((r).c[2],(s1).c[2],(s2).c[2]);\
 _vector_add_f((r).c[3],(s1).c[3],(s2).c[3])

#define _spinor_sub_f(r,s1,s2)\
 _vector_sub_f((r).c[0],(s1).c[0],(s2).c[0]);\
 _vector_sub_f((r).c[1],(s1).c[1],(s2).c[1]);\
 _vector_sub_f((r).c[2],(s1).c[2],(s2).c[2]);\
 _vector_sub_f((r).c[3],(s1).c[3],(s2).c[3])

#define _spinor_add_assign_f(r,s)\
 _vector_add_assign_f((r).c[0],(s).c[0]);\
 _vector_add_assign_f((r).c[1],(s).c[1]);\
 _vector_add_assign_f((r).c[2],(s).c[2]);\
 _vector_add_assign_f((r).c[3],(s).c[3])

#define _spinor_sub_assign_f(r,s)\
 _vector_sub_assign_f((r).c[0],(s).c[0]);\
 _vector_sub_assign_f((r).c[1],(s).c[1]);\
 _vector_sub_assign_f((r).c[2],(s).c[2]);\
 _vector_sub_assign_f((r).c[3],(s).c[3])

#define _spinor_prod_re_f(k,r,s)\
 do { _vector_prod_re_f((k),(r).c[0],(s).c[0]);\
 _vector_prod_add_assign_re_f((k),(r).c[1],(s).c[1]);\
 _vector_prod_add_assign_re_f((k),(r).c[2],(s).c[2]);\
 _vector_prod_add_assign_re_f((k),(r).c[3],(s).c[3]);\
 } while(0)

#define _spinor_prod_im_f(k,r,s)\
 do { _vector_prod_im_f((k),(r).c[0],(s).c[0]);\
 _vector_prod_add_assign_im_f((k),(r).c[1],(s).c[1]);\
 _vector_prod_add_assign_im_f((k),(r).c[2],(s).c[2]);\
 _vector_prod_add_assign_im_f((k),(r).c[3],(s).c[3]);\
 } while(0)

#define _spinor_prod_f(z,r,s)\
 do { _complex_0(z);\
 _vector_prod_assign_f((z),(r).c[0],(s).c[0]);\
 _vector_prod_assign_f((z),(r).c[1],(s).c[1]);\
 _vector_prod_assign_f((z),(r).c[2],(s).c[2]);\
 _vector_prod_assign_f((z),(r).c[3],(s).c[3]);\
 } while(0)

#define _spinor_g5_prod_re_f(k,r,s)\
 do { _vector_prod_re_f((k),(r).c[0],(s).c[0]);\
 _vector_prod_add_assign_re_f((k),(r).c[1],(s).c[1]);\
 _vector_prod_sub_assign_re_f((k),(r).c[2],(s).c[2]);\
 _vector_prod_sub_assign_re_f((k),(r).c[3],(s).c[3]);\
 } while(0)

#define _spinor_g5_prod_im_f(k,r,s)\
 do { _vector_prod_im_f((k),(r).c[0],(s).c[0]);\
 _vector_prod_add_assign_im_f((k),(r).c[1],(s).c[1]);\
 _vector_prod_sub_assign_im_f((k),(r).c[2],(s).c[2]);\
 _vector_prod_sub_assign_im_f((k),(r).c[3],(s).c[3]);\
 } while(0)

