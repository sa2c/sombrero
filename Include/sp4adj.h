#define _vector_zero_g(r)\
 _complex_0((r).c[0]);\
 _complex_0((r).c[1]);\
 _complex_0((r).c[2]);\
 _complex_0((r).c[3])

#define _algebra_vector_zero_g(r)\
 (r).c[0]=0.;\
 (r).c[1]=0.;\
 (r).c[2]=0.;\
 (r).c[3]=0.;\
 (r).c[4]=0.;\
 (r).c[5]=0.;\
 (r).c[6]=0.;\
 (r).c[7]=0.;\
 (r).c[8]=0.;\
 (r).c[9]=0.

#define _suNg_zero(u)\
 _complex_0((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3]);\
 _complex_0((u).c[4]);\
 _complex_0((u).c[5]);\
 _complex_0((u).c[6]);\
 _complex_0((u).c[7])

#define _suNg_unit(u)\
 _complex_1((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3]);\
 _complex_0((u).c[4]);\
 _complex_1((u).c[5]);\
 _complex_0((u).c[6]);\
 _complex_0((u).c[7]);\


#define _suNg_expand(v,u)\
 (v).c[0] = _complex_re((u).c[0]);\
 (v).c[0]+= I*( _complex_im((u).c[0]));\
 (v).c[10] = _complex_re((u).c[0]);\
 (v).c[10]+= I*( -_complex_im((u).c[0]));\
 (v).c[2] = _complex_re((u).c[2]);\
 (v).c[2]+= I*( _complex_im((u).c[2]));\
 (v).c[8] = -_complex_re((u).c[2]);\
 (v).c[8]+= I*( _complex_im((u).c[2]));\
 (v).c[1] = _complex_re((u).c[1]);\
 (v).c[1]+= I*( _complex_im((u).c[1]));\
 (v).c[11] = _complex_re((u).c[1]);\
 (v).c[11]+= I*( -_complex_im((u).c[1]));\
 (v).c[3] = _complex_re((u).c[3]);\
 (v).c[3]+= I*( _complex_im((u).c[3]));\
 (v).c[9] = -_complex_re((u).c[3]);\
 (v).c[9]+= I*( _complex_im((u).c[3]));\
 (v).c[4] = _complex_re((u).c[4]);\
 (v).c[4]+= I*( _complex_im((u).c[4]));\
 (v).c[14] = _complex_re((u).c[4]);\
 (v).c[14]+= I*( -_complex_im((u).c[4]));\
 (v).c[6] = _complex_re((u).c[6]);\
 (v).c[6]+= I*( _complex_im((u).c[6]));\
 (v).c[12] = -_complex_re((u).c[6]);\
 (v).c[12]+= I*( _complex_im((u).c[6]));\
 (v).c[5] = _complex_re((u).c[5]);\
 (v).c[5]+= I*( _complex_im((u).c[5]));\
 (v).c[15] = _complex_re((u).c[5]);\
 (v).c[15]+= I*( -_complex_im((u).c[5]));\
 (v).c[7] = _complex_re((u).c[7]);\
 (v).c[7]+= I*( _complex_im((u).c[7]));\
 (v).c[13] = -_complex_re((u).c[7]);\
 (v).c[13]+= I*( _complex_im((u).c[7]));\


#define _suNgfull_unit(u)\
 _complex_1((u).c[0]);\
 _complex_0((u).c[1]);\
 _complex_0((u).c[2]);\
 _complex_0((u).c[3]);\
 _complex_0((u).c[4]);\
 _complex_1((u).c[5]);\
 _complex_0((u).c[6]);\
 _complex_0((u).c[7]);\
 _complex_0((u).c[8]);\
 _complex_0((u).c[9]);\
 _complex_1((u).c[10]);\
 _complex_0((u).c[11]);\
 _complex_0((u).c[12]);\
 _complex_0((u).c[13]);\
 _complex_0((u).c[14]);\
 _complex_1((u).c[15])

#define _vector_zero_f(r)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_0((r).c[_i]);\
 ++_i;\
 _complex_0((r).c[_i]);\
 ++_i;\
 _complex_0((r).c[_i]);\
 ++_i;\
 _complex_0((r).c[_i]);\
 ++_i;\
 } _complex_0((r).c[_i]);\
 ++_i;\
 _complex_0((r).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_minus_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_plus_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_plus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_minus_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_minus((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_mul_f(r,k,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 } _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_mulc_f(r,z,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 } _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_mulc_star_f(r,z,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 } _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 _complex_mul_star((r).c[_i],(s).c[_i],(z));\
 ++_i;\
 } while(0)

#define _vector_add_f(r,s1,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_sub_f(r,s1,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_add_f(r,s1,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_sub_f(r,s1,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_add_assign_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_sub_assign_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_add_assign_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_add_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_i_sub_assign_f(r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_i_sub_assign((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_re_f(k,r,s)\
 do { int _i;\
 (k)=_complex_prod_re((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_im_f(k,r,s)\
 do { int _i;\
 (k)=_complex_prod_im((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_mulc_add_assign_f(r,z,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 } _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 _complex_mul_assign((r).c[_i],(z),(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_mul_add_assign_f(r,k,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 } _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 _complex_mulr_assign((r).c[_i],(k),(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_lc_f(r,k1,s1,k2,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 } _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_lc_add_assign_f(r,k1,s1,k2,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 } _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_clc_f(r,z1,s1,z2,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 } _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_clc_add_assign_f(r,z1,s1,z2,s2)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 } _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_assign_f(z,r,s)\
 do { int _i;\
for (_i=0;\
 _i<8;\
 ){ _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 } _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 _complex_prod_assign((z),(r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_add_assign_re_f(k,r,s)\
 do { int _i;\
 (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)+=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_add_assign_im_f(k,r,s)\
 do { int _i;\
 (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)+=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_sub_assign_re_f(k,r,s)\
 do { int _i;\
 (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)-=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)-=_complex_prod_re((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _vector_prod_sub_assign_im_f(k,r,s)\
 do { int _i;\
 (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
 for (_i=1;\
 _i<8;\
 ){ (k)-=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 (k)-=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } (k)-=_complex_prod_im((r).c[_i],(s).c[_i]);\
 ++_i;\
 } while(0)

#define _suNfc_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<10;\
 ++_i){ _complex_mul((r).c[_i],(u).c[_k],(s).c[0]);\
 ++_k;\
 int _j=1;\
 for (;\
 _j<8;\
 ){ _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 } _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 } } while(0)

#define _suNfc_inverse_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<10;\
 ++_i){ _complex_mul_star((r).c[_i],(s).c[0],(u).c[_k]);\
 int _j=1;\
 for (;\
 _j<8;\
 ){ _k+=10;\
 _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]);\
 ++_j;\
 _k+=10;\
 _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]);\
 ++_j;\
 _k+=10;\
 _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]);\
 ++_j;\
 _k+=10;\
 _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]);\
 ++_j;\
 } _k+=10;\
 _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]);\
 ++_j;\
 _k-=89;\
 } } while(0)

#define _suNfc_zero(u)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 _complex_0((u).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNf_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<10;\
 ++_i){ _complex_mulr((r).c[_i],(u).c[_k],(s).c[0]);\
 ++_k;\
 int _j=1;\
 for (;\
 _j<8;\
 ){ _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 } _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_k;\
 ++_j;\
 } } while(0)

#define _suNf_inverse_multiply(r,u,s)\
 do { int _i,_k=0;\
for (_i=0;\
 _i<10;\
 ++_i){ _complex_mulr((r).c[_i],(u).c[_k],(s).c[0]);\
 int _j=1;\
 for (;\
 _j<8;\
 ){ _k+=10;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_j;\
 _k+=10;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_j;\
 _k+=10;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_j;\
 _k+=10;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_j;\
 } _k+=10;\
 _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]);\
 ++_j;\
 _k-=89;\
 } } while(0)

#define _suNf_zero(u)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]=0.;\
 ++_i;\
 (u).c[_i]=0.;\
 ++_i;\
 (u).c[_i]=0.;\
 ++_i;\
 (u).c[_i]=0.;\
 ++_i;\
 } } while(0)

#define _suNfc_dagger(u,v)\
 do { int _i,_j,_n=0,_k=0;\
 for (_i=0;\
 _i<10;\
 ++_i){ _complex_star((u).c[_n],(v).c[_k]);\
 for (_j=0;\
 _j<8;\
 ){ ++_n;\
 _k+=10;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_j;\
 ++_n;\
 _k+=10;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_j;\
 ++_n;\
 _k+=10;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_j;\
 ++_n;\
 _k+=10;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_j;\
 } ++_n;\
 _k+=10;\
 _complex_star((u).c[_n],(v).c[_k]);\
 ++_n;\
 _k-=89;\
 } } while(0)

#define _suNfc_mul(u,r,v)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 _complex_mulr((u).c[_i],(r),(v).c[_i]);\
 ++_i;\
 } } while(0)

#define _suNfc_mul_assign(u,r)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 } } while(0)

#define _suNf_dagger(u,v)\
 { int _i,_j,_n=0,_k=0;\
 for (_i=0;\
 _i<10;\
 ++_i){ (u).c[_n]=(v).c[_k];\
 for (_j=0;\
 _j<8;\
 ){ ++_n;\
 _k+=10;\
 (u).c[_n]=(v).c[_k];\
 ++_j;\
 ++_n;\
 _k+=10;\
 (u).c[_n]=(v).c[_k];\
 ++_j;\
 ++_n;\
 _k+=10;\
 (u).c[_n]=(v).c[_k];\
 ++_j;\
 ++_n;\
 _k+=10;\
 (u).c[_n]=(v).c[_k];\
 ++_j;\
 } ++_n;\
 _k+=10;\
 (u).c[_n]=(v).c[_k];\
 ++_n;\
 _k-=89;\
 } }((void)0)

#define _suNf_times_suNf(u,v,w)\
 do { int _i,_y,_j,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<10;\
 ++_i){ for (_y=0;\
 _y<10;\
 ++_y){ (u).c[_n]=(v).c[_k]*(w).c[_l];\
 for (_j=0;\
 _j<8;\
 ){ ++_k;\
 _l+=10;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 _l+=10;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 _l+=10;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 _l+=10;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 } ++_k;\
 _l+=10;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_n;\
 _k-=9;\
 _l-=89;\
 } _k+=10;\
 _l=0;\
 } } while(0)

#define _suNf_times_suNf_dagger(u,v,w)\
 do { int _i,_y,_j,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<10;\
 ++_i){ for (_y=0;\
 _y<10;\
 ++_y){ (u).c[_n]=(v).c[_k]*(w).c[_l];\
 for (_j=0;\
 _j<8;\
 ){ ++_k;\
 ++_l;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 ++_l;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 ++_l;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 ++_k;\
 ++_l;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_j;\
 } ++_k;\
 ++_l;\
 (u).c[_n]+=(v).c[_k]*(w).c[_l];\
 ++_n;\
 _k-=9;\
 ++_l;\
 } _k+=10;\
 _l=0;\
 } } while(0)

#define _suNf_dagger_times_suNf(u,v,w)\
 do { int _i,_y,_j,_n=0,_k=0,_l=0;\
 for (_i=0;\
 _i<10;\
 ++_i){ for (_y=0;\
 _y<10;\
 ++_y){ _k=_y;\
 _l=_i;\
 (u).c[_n]=(w).c[_k]*(v).c[_l];\
 for (_j=0;\
 _j<8;\
 ){ _k+=10;\
 _l+=10;\
 (u).c[_n]+=(w).c[_k]*(v).c[_l];\
 ++_j;\
 _k+=10;\
 _l+=10;\
 (u).c[_n]+=(w).c[_k]*(v).c[_l];\
 ++_j;\
 _k+=10;\
 _l+=10;\
 (u).c[_n]+=(w).c[_k]*(v).c[_l];\
 ++_j;\
 _k+=10;\
 _l+=10;\
 (u).c[_n]+=(w).c[_k]*(v).c[_l];\
 ++_j;\
 } _k+=10;\
 _l+=10;\
 (u).c[_n]+=(w).c[_k]*(v).c[_l];\
 ++_n;\
 } } } while(0)

#define _suNf_add_assign(u,v)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]+=(v).c[_i];\
 ++_i;\
 (u).c[_i]+=(v).c[_i];\
 ++_i;\
 (u).c[_i]+=(v).c[_i];\
 ++_i;\
 (u).c[_i]+=(v).c[_i];\
 ++_i;\
 } } while(0)

#define _suNf_mul(u,r,v)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]=(r)*(v).c[_i];\
 ++_i;\
 (u).c[_i]=(r)*(v).c[_i];\
 ++_i;\
 (u).c[_i]=(r)*(v).c[_i];\
 ++_i;\
 (u).c[_i]=(r)*(v).c[_i];\
 ++_i;\
 } } while(0)

#define _suNf_mul_assign(u,r)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 (u).c[_i]*=(r);\
 ++_i;\
 } } while(0)

#define _suNf_minus(u,v)\
 do { int _i;\
for (_i=0;\
 _i<100;\
 ){ (u).c[_i]=-(v).c[_i];\
 ++_i;\
 (u).c[_i]=-(v).c[_i];\
 ++_i;\
 (u).c[_i]=-(v).c[_i];\
 ++_i;\
 (u).c[_i]=-(v).c[_i];\
 ++_i;\
 } } while(0)

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

