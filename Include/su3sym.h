#define _vector_zero_g(r)\
 _complex_0((r).c[0]);\
 _complex_0((r).c[1]);\
 _complex_0((r).c[2])

#define _suNg_zero(u)\
 _complex_0((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3]);\
 _complex_0((u).c[4]);\
 _complex_0((u).c[5]);\
 _complex_0((u).c[6]);\
 _complex_0((u).c[7]);\
 _complex_0((u).c[8])

#define _algebra_vector_zero_g(r)\
 (r).c[0]=0.;\
 (r).c[1]=0.;\
 (r).c[2]=0.;\
 (r).c[3]=0.;\
 (r).c[4]=0.;\
 (r).c[5]=0.;\
 (r).c[6]=0.;\
 (r).c[7]=0.

#define _suNg_unit(u)\
 _complex_1((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3]);\
 _complex_1((u).c[4]);\
 _complex_0((u).c[5]);\
 _complex_0((u).c[6]);\
 _complex_0((u).c[7]);\
 _complex_1((u).c[8])

#define _vector_zero_f(r)\
 _complex_0((r).c[0]);\
 _complex_0((r).c[1]);\
 _complex_0((r).c[2]);\
 _complex_0((r).c[3]);\
 _complex_0((r).c[4]);\
 _complex_0((r).c[5])

#define _vector_minus_f(r,s)\
 _complex_minus((r).c[0],(s).c[0]);\
 _complex_minus((r).c[1],(s).c[1]);\
 _complex_minus((r).c[2],(s).c[2]);\
 _complex_minus((r).c[3],(s).c[3]);\
 _complex_minus((r).c[4],(s).c[4]);\
 _complex_minus((r).c[5],(s).c[5])

#define _vector_i_plus_f(r,s)\
 _complex_i_plus((r).c[0],(s).c[0]);\
 _complex_i_plus((r).c[1],(s).c[1]);\
 _complex_i_plus((r).c[2],(s).c[2]);\
 _complex_i_plus((r).c[3],(s).c[3]);\
 _complex_i_plus((r).c[4],(s).c[4]);\
 _complex_i_plus((r).c[5],(s).c[5])

#define _vector_i_minus_f(r,s)\
 _complex_i_minus((r).c[0],(s).c[0]);\
 _complex_i_minus((r).c[1],(s).c[1]);\
 _complex_i_minus((r).c[2],(s).c[2]);\
 _complex_i_minus((r).c[3],(s).c[3]);\
 _complex_i_minus((r).c[4],(s).c[4]);\
 _complex_i_minus((r).c[5],(s).c[5])

#define _vector_mul_f(r,k,s)\
 _complex_mulr((r).c[0],(k),(s).c[0]);\
 _complex_mulr((r).c[1],(k),(s).c[1]);\
 _complex_mulr((r).c[2],(k),(s).c[2]);\
 _complex_mulr((r).c[3],(k),(s).c[3]);\
 _complex_mulr((r).c[4],(k),(s).c[4]);\
 _complex_mulr((r).c[5],(k),(s).c[5])

#define _vector_mulc_f(r,z,s)\
 _complex_mul((r).c[0],(z),(s).c[0]);\
 _complex_mul((r).c[1],(z),(s).c[1]);\
 _complex_mul((r).c[2],(z),(s).c[2]);\
 _complex_mul((r).c[3],(z),(s).c[3]);\
 _complex_mul((r).c[4],(z),(s).c[4]);\
 _complex_mul((r).c[5],(z),(s).c[5])

#define _vector_mulc_star_f(r,z,s)\
 _complex_mul_star((r).c[0],(s).c[0],(z));\
 _complex_mul_star((r).c[1],(s).c[1],(z));\
 _complex_mul_star((r).c[2],(s).c[2],(z));\
 _complex_mul_star((r).c[3],(s).c[3],(z));\
 _complex_mul_star((r).c[4],(s).c[4],(z));\
 _complex_mul_star((r).c[5],(s).c[5],(z))

#define _vector_add_f(r,s1,s2)\
 _complex_add((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_add((r).c[1],(s1).c[1],(s2).c[1]);\
 _complex_add((r).c[2],(s1).c[2],(s2).c[2]);\
 _complex_add((r).c[3],(s1).c[3],(s2).c[3]);\
 _complex_add((r).c[4],(s1).c[4],(s2).c[4]);\
 _complex_add((r).c[5],(s1).c[5],(s2).c[5])

#define _vector_sub_f(r,s1,s2)\
 _complex_sub((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_sub((r).c[1],(s1).c[1],(s2).c[1]);\
 _complex_sub((r).c[2],(s1).c[2],(s2).c[2]);\
 _complex_sub((r).c[3],(s1).c[3],(s2).c[3]);\
 _complex_sub((r).c[4],(s1).c[4],(s2).c[4]);\
 _complex_sub((r).c[5],(s1).c[5],(s2).c[5])

#define _vector_i_add_f(r,s1,s2)\
 _complex_i_add((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_i_add((r).c[1],(s1).c[1],(s2).c[1]);\
 _complex_i_add((r).c[2],(s1).c[2],(s2).c[2]);\
 _complex_i_add((r).c[3],(s1).c[3],(s2).c[3]);\
 _complex_i_add((r).c[4],(s1).c[4],(s2).c[4]);\
 _complex_i_add((r).c[5],(s1).c[5],(s2).c[5])

#define _vector_i_sub_f(r,s1,s2)\
 _complex_i_sub((r).c[0],(s1).c[0],(s2).c[0]);\
 _complex_i_sub((r).c[1],(s1).c[1],(s2).c[1]);\
 _complex_i_sub((r).c[2],(s1).c[2],(s2).c[2]);\
 _complex_i_sub((r).c[3],(s1).c[3],(s2).c[3]);\
 _complex_i_sub((r).c[4],(s1).c[4],(s2).c[4]);\
 _complex_i_sub((r).c[5],(s1).c[5],(s2).c[5])

#define _vector_add_assign_f(r,s)\
 _complex_add_assign((r).c[0],(s).c[0]);\
 _complex_add_assign((r).c[1],(s).c[1]);\
 _complex_add_assign((r).c[2],(s).c[2]);\
 _complex_add_assign((r).c[3],(s).c[3]);\
 _complex_add_assign((r).c[4],(s).c[4]);\
 _complex_add_assign((r).c[5],(s).c[5])

#define _vector_sub_assign_f(r,s)\
 _complex_sub_assign((r).c[0],(s).c[0]);\
 _complex_sub_assign((r).c[1],(s).c[1]);\
 _complex_sub_assign((r).c[2],(s).c[2]);\
 _complex_sub_assign((r).c[3],(s).c[3]);\
 _complex_sub_assign((r).c[4],(s).c[4]);\
 _complex_sub_assign((r).c[5],(s).c[5])

#define _vector_i_add_assign_f(r,s)\
 _complex_i_add_assign((r).c[0],(s).c[0]);\
 _complex_i_add_assign((r).c[1],(s).c[1]);\
 _complex_i_add_assign((r).c[2],(s).c[2]);\
 _complex_i_add_assign((r).c[3],(s).c[3]);\
 _complex_i_add_assign((r).c[4],(s).c[4]);\
 _complex_i_add_assign((r).c[5],(s).c[5])

#define _vector_i_sub_assign_f(r,s)\
 _complex_i_sub_assign((r).c[0],(s).c[0]);\
 _complex_i_sub_assign((r).c[1],(s).c[1]);\
 _complex_i_sub_assign((r).c[2],(s).c[2]);\
 _complex_i_sub_assign((r).c[3],(s).c[3]);\
 _complex_i_sub_assign((r).c[4],(s).c[4]);\
 _complex_i_sub_assign((r).c[5],(s).c[5])

#define _vector_prod_re_f(k,r,s)\
 (k)=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_re((r).c[1],(s).c[1]);\
 (k)+=_complex_prod_re((r).c[2],(s).c[2]);\
 (k)+=_complex_prod_re((r).c[3],(s).c[3]);\
 (k)+=_complex_prod_re((r).c[4],(s).c[4]);\
 (k)+=_complex_prod_re((r).c[5],(s).c[5])

#define _vector_prod_im_f(k,r,s)\
 (k)=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_im((r).c[1],(s).c[1]);\
 (k)+=_complex_prod_im((r).c[2],(s).c[2]);\
 (k)+=_complex_prod_im((r).c[3],(s).c[3]);\
 (k)+=_complex_prod_im((r).c[4],(s).c[4]);\
 (k)+=_complex_prod_im((r).c[5],(s).c[5])

#define _vector_mulc_add_assign_f(r,z,s)\
 _complex_mul_assign((r).c[0],(z),(s).c[0]);\
 _complex_mul_assign((r).c[1],(z),(s).c[1]);\
 _complex_mul_assign((r).c[2],(z),(s).c[2]);\
 _complex_mul_assign((r).c[3],(z),(s).c[3]);\
 _complex_mul_assign((r).c[4],(z),(s).c[4]);\
 _complex_mul_assign((r).c[5],(z),(s).c[5])

#define _vector_mul_add_assign_f(r,k,s)\
 _complex_mulr_assign((r).c[0],(k),(s).c[0]);\
 _complex_mulr_assign((r).c[1],(k),(s).c[1]);\
 _complex_mulr_assign((r).c[2],(k),(s).c[2]);\
 _complex_mulr_assign((r).c[3],(k),(s).c[3]);\
 _complex_mulr_assign((r).c[4],(k),(s).c[4]);\
 _complex_mulr_assign((r).c[5],(k),(s).c[5])

#define _vector_lc_f(r,k1,s1,k2,s2)\
 _complex_rlc((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _complex_rlc((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]);\
 _complex_rlc((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]);\
 _complex_rlc((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3]);\
 _complex_rlc((r).c[4],(k1),(s1).c[4],(k2),(s2).c[4]);\
 _complex_rlc((r).c[5],(k1),(s1).c[5],(k2),(s2).c[5])

#define _vector_lc_add_assign_f(r,k1,s1,k2,s2)\
 _complex_rlc_assign((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]);\
 _complex_rlc_assign((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]);\
 _complex_rlc_assign((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]);\
 _complex_rlc_assign((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3]);\
 _complex_rlc_assign((r).c[4],(k1),(s1).c[4],(k2),(s2).c[4]);\
 _complex_rlc_assign((r).c[5],(k1),(s1).c[5],(k2),(s2).c[5])

#define _vector_clc_f(r,z1,s1,z2,s2)\
 _complex_clc((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _complex_clc((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]);\
 _complex_clc((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]);\
 _complex_clc((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3]);\
 _complex_clc((r).c[4],(z1),(s1).c[4],(z2),(s2).c[4]);\
 _complex_clc((r).c[5],(z1),(s1).c[5],(z2),(s2).c[5])

#define _vector_clc_add_assign_f(r,z1,s1,z2,s2)\
 _complex_clc_assign((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]);\
 _complex_clc_assign((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]);\
 _complex_clc_assign((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]);\
 _complex_clc_assign((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3]);\
 _complex_clc_assign((r).c[4],(z1),(s1).c[4],(z2),(s2).c[4]);\
 _complex_clc_assign((r).c[5],(z1),(s1).c[5],(z2),(s2).c[5])

#define _vector_prod_assign_f(z,r,s)\
 _complex_prod_assign((z),(r).c[0],(s).c[0]);\
 _complex_prod_assign((z),(r).c[1],(s).c[1]);\
 _complex_prod_assign((z),(r).c[2],(s).c[2]);\
 _complex_prod_assign((z),(r).c[3],(s).c[3]);\
 _complex_prod_assign((z),(r).c[4],(s).c[4]);\
 _complex_prod_assign((z),(r).c[5],(s).c[5])

#define _vector_prod_add_assign_re_f(k,r,s)\
 (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_re((r).c[1],(s).c[1]);\
 (k)+=_complex_prod_re((r).c[2],(s).c[2]);\
 (k)+=_complex_prod_re((r).c[3],(s).c[3]);\
 (k)+=_complex_prod_re((r).c[4],(s).c[4]);\
 (k)+=_complex_prod_re((r).c[5],(s).c[5])

#define _vector_prod_add_assign_im_f(k,r,s)\
 (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)+=_complex_prod_im((r).c[1],(s).c[1]);\
 (k)+=_complex_prod_im((r).c[2],(s).c[2]);\
 (k)+=_complex_prod_im((r).c[3],(s).c[3]);\
 (k)+=_complex_prod_im((r).c[4],(s).c[4]);\
 (k)+=_complex_prod_im((r).c[5],(s).c[5])

#define _vector_prod_sub_assign_re_f(k,r,s)\
 (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
 (k)-=_complex_prod_re((r).c[1],(s).c[1]);\
 (k)-=_complex_prod_re((r).c[2],(s).c[2]);\
 (k)-=_complex_prod_re((r).c[3],(s).c[3]);\
 (k)-=_complex_prod_re((r).c[4],(s).c[4]);\
 (k)-=_complex_prod_re((r).c[5],(s).c[5])

#define _vector_prod_sub_assign_im_f(k,r,s)\
 (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
 (k)-=_complex_prod_im((r).c[1],(s).c[1]);\
 (k)-=_complex_prod_im((r).c[2],(s).c[2]);\
 (k)-=_complex_prod_im((r).c[3],(s).c[3]);\
 (k)-=_complex_prod_im((r).c[4],(s).c[4]);\
 (k)-=_complex_prod_im((r).c[5],(s).c[5])

#define _suNf_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<6;\
 ++_i){ _complex_mul((r).c[_i],(u).c[_k],(s).c[0]);\
 ++_k;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[1]);\
 ++_k;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[2]);\
 ++_k;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[3]);\
 ++_k;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[4]);\
 ++_k;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[5]);\
 ++_k;\
 } } while(0)

#define _suNf_inverse_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<6;\
 ++_i){ _complex_mul_star((r).c[_i],(s).c[0],(u).c[_k]);\
 _k+=6;\
 _complex_mul_star_assign((r).c[_i],(s).c[1],(u).c[_k]);\
 _k+=6;\
 _complex_mul_star_assign((r).c[_i],(s).c[2],(u).c[_k]);\
 _k+=6;\
 _complex_mul_star_assign((r).c[_i],(s).c[3],(u).c[_k]);\
 _k+=6;\
 _complex_mul_star_assign((r).c[_i],(s).c[4],(u).c[_k]);\
 _k+=6;\
 _complex_mul_star_assign((r).c[_i],(s).c[5],(u).c[_k]);\
 _k-=29;\
 } } while(0)

#define _suNf_zero(u)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNf_dagger(u,v)\
 do { int _i,_n=0,_k=0;\
 for (_i=0;\
 _i<6;\
 ++_i){ _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k+=6;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k+=6;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k+=6;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k+=6;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k+=6;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k-=29;\
 } } while(0)

#define _suNf_times_suNf(u,v,w)\
 do { int _i,_y,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<6;\
 ++_i){ for (_y=0;\
 _y<6;\
 ++_y){ _complex_mul((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 _l+=6;\
 _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 _l+=6;\
 _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 _l+=6;\
 _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 _l+=6;\
 _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 _l+=6;\
 _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_n;\
 _k-=5;\
 _l-=29;\
 } _k+=6;\
 _l=0;\
 } } while(0)

#define _suNf_times_suNf_dagger(u,v,w)\
 do { int _i,_y,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<6;\
 ++_i){ for (_y=0;\
 _y<6;\
 ++_y){ _complex_mul_star((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 ++_l;\
 _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 ++_l;\
 _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 ++_l;\
 _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 ++_l;\
 _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_k;\
 ++_l;\
 _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
 ++_n;\
 _k-=5;\
 ++_l;\
 } _k+=6;\
 _l=0;\
 } } while(0)

#define _suNf_dagger_times_suNf(u,v,w)\
 do { int _i,_y,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<6;\
 ++_i){ for (_y=0;\
 _y<6;\
 ++_y){ _k=_y;\
 _l=_i;\
 _complex_mul_star((u).c[_n],(w).c[_k],(v).c[_l]);\
 _k+=6;\
 _l+=6;\
 _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
 _k+=6;\
 _l+=6;\
 _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
 _k+=6;\
 _l+=6;\
 _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
 _k+=6;\
 _l+=6;\
 _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
 _k+=6;\
 _l+=6;\
 _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
 ++_n;\
 } } } while(0)

#define _suNf_minus(u,v)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ _complex_minus((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_minus((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_minus((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_minus((u).c[_i],(v).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNf_mul(u,r,v)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNf_mul_assign(u,r)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 } } while(0)

#define _suNf_mulc(u,r,v)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ _complex_mul((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mul((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mul((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mul((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNf_add_assign(u,v)\
 do { int _i;\
for (_i=0;\
 _i<36;\
 ){ _complex_add_assign((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_add_assign((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_add_assign((u).c[_i],(v).c[_i]);\
 ++_i;\
 _complex_add_assign((u).c[_i],(v).c[_i]);\
 ++_i;\
 } } while(0)

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

