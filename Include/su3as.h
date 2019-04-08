/*******************************************************************************
*
* File suN.h
*
* Type definitions and macros for SU(N) matrices and spinors
*
*******************************************************************************/

#ifndef SUN_H
#define SUN_H

#include "suN_types.h"

/*******************************************************************************
*
* The following macros are the same for single and double precision types
*
* Depending on the macro, arguments are variables of type suN_vector and suN
* (or suN_vector_flt and suN_flt)
*
*******************************************************************************/

/* r=0 */
#define _vector_zero_g(r) \
   _complex_0((r).c[0]); \
   _complex_0((r).c[1]); \
   _complex_0((r).c[2])

/* r=-s */
#define _vector_minus_g(r,s) \
   _complex_minus((r).c[0],(s).c[0]); \
   _complex_minus((r).c[1],(s).c[1]); \
   _complex_minus((r).c[2],(s).c[2])

/* r= i*s */
#define _vector_i_plus_g(r,s) \
   _complex_i_plus((r).c[0],(s).c[0]); \
   _complex_i_plus((r).c[1],(s).c[1]); \
   _complex_i_plus((r).c[2],(s).c[2])

/* r=-i*s */
#define _vector_i_minus_g(r,s) \
   _complex_i_minus((r).c[0],(s).c[0]); \
   _complex_i_minus((r).c[1],(s).c[1]); \
   _complex_i_minus((r).c[2],(s).c[2])

/* r=k*s (k real) */
#define _vector_mul_g(r,k,s) \
   _complex_mulr((r).c[0],(k),(s).c[0]); \
   _complex_mulr((r).c[1],(k),(s).c[1]); \
   _complex_mulr((r).c[2],(k),(s).c[2])

/* r=z*s (z complex) */
#define _vector_mulc_g(r,z,s) \
   _complex_mul((r).c[0],(z),(s).c[0]); \
   _complex_mul((r).c[1],(z),(s).c[1]); \
   _complex_mul((r).c[2],(z),(s).c[2])

/* r=(z^+)*s (z complex) */
#define _vector_mulc_star_g(r,z,s) \
   _complex_mul_star((r).c[0],(s).c[0],(z)); \
   _complex_mul_star((r).c[1],(s).c[1],(z)); \
   _complex_mul_star((r).c[2],(s).c[2],(z))

/* r=s1+s2 */
#define _vector_add_g(r,s1,s2) \
   _complex_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_add((r).c[2],(s1).c[2],(s2).c[2])

/* r=s1-s2 */
#define _vector_sub_g(r,s1,s2) \
   _complex_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_sub((r).c[2],(s1).c[2],(s2).c[2])

/* r=s1+i*s2 */
#define _vector_i_add_g(r,s1,s2) \
   _complex_i_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_add((r).c[2],(s1).c[2],(s2).c[2])

/* r=s1-i*s2 */
#define _vector_i_sub_g(r,s1,s2) \
   _complex_i_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_sub((r).c[2],(s1).c[2],(s2).c[2])

/* r+=s */
#define _vector_add_assign_g(r,s) \
   _complex_add_assign((r).c[0],(s).c[0]); \
   _complex_add_assign((r).c[1],(s).c[1]); \
   _complex_add_assign((r).c[2],(s).c[2])

/* r-=s */
#define _vector_sub_assign_g(r,s) \
   _complex_sub_assign((r).c[0],(s).c[0]); \
   _complex_sub_assign((r).c[1],(s).c[1]); \
   _complex_sub_assign((r).c[2],(s).c[2])

/* r+=i*s */
#define _vector_i_add_assign_g(r,s) \
   _complex_i_add_assign((r).c[0],(s).c[0]); \
   _complex_i_add_assign((r).c[1],(s).c[1]); \
   _complex_i_add_assign((r).c[2],(s).c[2])

/* r-=i*s */
#define _vector_i_sub_assign_g(r,s) \
   _complex_i_sub_assign((r).c[0],(s).c[0]); \
   _complex_i_sub_assign((r).c[1],(s).c[1]); \
   _complex_i_sub_assign((r).c[2],(s).c[2])

/* k=Re(r^*s) */
#define _vector_prod_re_g(k,r,s) \
   (k)=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2])

/* k=Im(r*s) */
#define _vector_prod_im_g(k,r,s) \
   (k)=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2])

/* r+=z*s (z complex) */
#define _vector_mulc_add_assign_g(r,z,s) \
   _complex_mul_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_assign((r).c[2],(z),(s).c[2])

/* r+=k*s (k real) */
#define _vector_mul_add_assign_g(r,k,s) \
   _complex_mulr_assign((r).c[0],(k),(s).c[0]); \
   _complex_mulr_assign((r).c[1],(k),(s).c[1]); \
   _complex_mulr_assign((r).c[2],(k),(s).c[2])

/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_g(r,k1,s1,k2,s2) \
   _complex_rlc((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2])

/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_add_assign_g(r,k1,s1,k2,s2) \
   _complex_rlc_assign((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc_assign((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc_assign((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_g(r,z1,s1,z2,s2) \
   _complex_clc((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_add_assign_g(r,z1,s1,z2,s2) \
   _complex_clc_assign((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc_assign((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc_assign((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2])

/* z+=r^*s (c complex) */
#define _vector_prod_assign_g(z,r,s) \
   _complex_prod_assign((z),(r).c[0],(s).c[0]); \
   _complex_prod_assign((z),(r).c[1],(s).c[1]); \
   _complex_prod_assign((z),(r).c[2],(s).c[2])

/* k+=Re(r^*s) */
#define _vector_prod_add_assign_re_g(k,r,s) \
   (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2])

/* k+=Im(r*s) */
#define _vector_prod_add_assign_im_g(k,r,s) \
   (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2])

/* k-=Re(r^*s) */
#define _vector_prod_sub_assign_re_g(k,r,s) \
   (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_re((r).c[2],(s).c[2])

/* k-=Im(r*s) */
#define _vector_prod_sub_assign_im_g(k,r,s) \
   (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_im((r).c[2],(s).c[2])

/* r-=z*s (z complex) */
#define _vector_project_g(r,z,s) \
   _complex_mul_sub_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_sub_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_sub_assign((r).c[2],(z),(s).c[2])

/* SU(N) matrix u times SU(N) vector s */
/* r=u*s */
#define _suNg_multiply(r,u,s) \
      _complex_mul((r).c[0],(u).c[0],(s).c[0]);\
      _complex_mul_assign((r).c[0],(u).c[1],(s).c[1]); \
      _complex_mul_assign((r).c[0],(u).c[2],(s).c[2]); \
      _complex_mul((r).c[1],(u).c[3],(s).c[0]);\
      _complex_mul_assign((r).c[1],(u).c[4],(s).c[1]); \
      _complex_mul_assign((r).c[1],(u).c[5],(s).c[2]); \
      _complex_mul((r).c[2],(u).c[6],(s).c[0]);\
      _complex_mul_assign((r).c[2],(u).c[7],(s).c[1]); \
      _complex_mul_assign((r).c[2],(u).c[8],(s).c[2])

/* SU(N) matrix u^dagger times SU(N) vector s */
/* r=(u^dagger)*s */
#define _suNg_inverse_multiply(r,u,s) \
      _complex_mul_star((r).c[0],(s).c[0],(u).c[0]);\
      _complex_mul_star_assign((r).c[0],(s).c[1],(u).c[3]); \
      _complex_mul_star_assign((r).c[0],(s).c[2],(u).c[6]); \
      _complex_mul_star((r).c[1],(s).c[0],(u).c[1]);\
      _complex_mul_star_assign((r).c[1],(s).c[1],(u).c[4]); \
      _complex_mul_star_assign((r).c[1],(s).c[2],(u).c[7]); \
      _complex_mul_star((r).c[2],(s).c[0],(u).c[2]);\
      _complex_mul_star_assign((r).c[2],(s).c[1],(u).c[5]); \
      _complex_mul_star_assign((r).c[2],(s).c[2],(u).c[8])

/* u=0 */
#define _suNg_zero(u) \
    _complex_0((u).c[0]);\
    _complex_0((u).c[1]);\
    _complex_0((u).c[2]);\
    _complex_0((u).c[3]);\
    _complex_0((u).c[4]);\
    _complex_0((u).c[5]);\
    _complex_0((u).c[6]);\
    _complex_0((u).c[7]);\
    _complex_0((u).c[8])

/* r+=s */
#define _algebra_vector_add_assign_g(r,s) \
      (r).c[0]+=(s).c[0]; \
      (r).c[1]+=(s).c[1]; \
      (r).c[2]+=(s).c[2]; \
      (r).c[3]+=(s).c[3]; \
      (r).c[4]+=(s).c[4]; \
      (r).c[5]+=(s).c[5]; \
      (r).c[6]+=(s).c[6]; \
      (r).c[7]+=(s).c[7]

/* r-=s */
#define _algebra_vector_sub_assign_g(r,s) \
      (r).c[0]-=(s).c[0]; \
      (r).c[1]-=(s).c[1]; \
      (r).c[2]-=(s).c[2]; \
      (r).c[3]-=(s).c[3]; \
      (r).c[4]-=(s).c[4]; \
      (r).c[5]-=(s).c[5]; \
      (r).c[6]-=(s).c[6]; \
      (r).c[7]-=(s).c[7]

/* r+=k*s (k real) */
#define _algebra_vector_mul_add_assign_g(r,k,s) \
      (r).c[0]+=(k)*(s).c[0]; \
      (r).c[1]+=(k)*(s).c[1]; \
      (r).c[2]+=(k)*(s).c[2]; \
      (r).c[3]+=(k)*(s).c[3]; \
      (r).c[4]+=(k)*(s).c[4]; \
      (r).c[5]+=(k)*(s).c[5]; \
      (r).c[6]+=(k)*(s).c[6]; \
      (r).c[7]+=(k)*(s).c[7]

/* r=k*s (k real) */
#define _algebra_vector_mul_g(r,k,s) \
      (r).c[0]=(k)*(s).c[0]; \
      (r).c[1]=(k)*(s).c[1]; \
      (r).c[2]=(k)*(s).c[2]; \
      (r).c[3]=(k)*(s).c[3]; \
      (r).c[4]=(k)*(s).c[4]; \
      (r).c[5]=(k)*(s).c[5]; \
      (r).c[6]=(k)*(s).c[6]; \
      (r).c[7]=(k)*(s).c[7]

/* r=0  */
#define _algebra_vector_zero_g(r) \
      (r).c[0]=0.; \
      (r).c[1]=0.; \
      (r).c[2]=0.; \
      (r).c[3]=0.; \
      (r).c[4]=0.; \
      (r).c[5]=0.; \
      (r).c[6]=0.; \
      (r).c[7]=0.

/* k=|v|^2  */
#define _algebra_vector_sqnorm_g(k,r) \
   (k)=((r).c[0]*(r).c[0])+ \
       ((r).c[1]*(r).c[1])+ \
       ((r).c[2]*(r).c[2])+ \
       ((r).c[3]*(r).c[3])+ \
       ((r).c[4]*(r).c[4])+ \
       ((r).c[5]*(r).c[5])+ \
       ((r).c[6]*(r).c[6])+ \
       ((r).c[7]*(r).c[7])

/* k=Scalar product r*s (r,s algabra vectors)  */
#define _algebra_vector_prod_g(k,r,s) \
   (k)=((r).c[0]*(s).c[0])+ \
       ((r).c[1]*(s).c[1])+ \
       ((r).c[2]*(s).c[2])+ \
       ((r).c[3]*(s).c[3])+ \
       ((r).c[4]*(s).c[4])+ \
       ((r).c[5]*(s).c[5])+ \
       ((r).c[6]*(s).c[6])+ \
       ((r).c[7]*(s).c[7])

/* u = Omega */
#define _symplectic(u) \
   _complex_0((u).c[0]);\
   _complex_0((u).c[1]);\
   _complex_0((u).c[2]);\
   _complex_0((u).c[3]);\
   _complex_0((u).c[4]);\
   _complex_0((u).c[5]);\
/* u = u* */
#define _vector_conjugate(u) \
    _complex_star_assign((u).c[0]);\
    _complex_star_assign((u).c[1]);\
    _complex_star_assign((u).c[2])

/*******************************************************************************
*
* Macros for SU(N) matrices
*
* Arguments are variables of type suN
*
*******************************************************************************/

/* u=v^dagger */
#define _suNg_dagger(u,v) \
   _complex_star((u).c[0],(v).c[0]); \
   _complex_star((u).c[1],(v).c[3]); \
   _complex_star((u).c[2],(v).c[6]); \
   _complex_star((u).c[3],(v).c[1]); \
   _complex_star((u).c[4],(v).c[4]); \
   _complex_star((u).c[5],(v).c[7]); \
   _complex_star((u).c[6],(v).c[2]); \
   _complex_star((u).c[7],(v).c[5]); \
   _complex_star((u).c[8],(v).c[8])

/* u=v*w */
#define _suNg_times_suNg(u,v,w) \
      _complex_mul((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_assign((u).c[0],(v).c[1],(w).c[3]); \
      _complex_mul_assign((u).c[0],(v).c[2],(w).c[6]); \
      _complex_mul((u).c[1],(v).c[0],(w).c[1]);\
      _complex_mul_assign((u).c[1],(v).c[1],(w).c[4]); \
      _complex_mul_assign((u).c[1],(v).c[2],(w).c[7]); \
      _complex_mul((u).c[2],(v).c[0],(w).c[2]);\
      _complex_mul_assign((u).c[2],(v).c[1],(w).c[5]); \
      _complex_mul_assign((u).c[2],(v).c[2],(w).c[8]); \
      _complex_mul((u).c[3],(v).c[3],(w).c[0]);\
      _complex_mul_assign((u).c[3],(v).c[4],(w).c[3]); \
      _complex_mul_assign((u).c[3],(v).c[5],(w).c[6]); \
      _complex_mul((u).c[4],(v).c[3],(w).c[1]);\
      _complex_mul_assign((u).c[4],(v).c[4],(w).c[4]); \
      _complex_mul_assign((u).c[4],(v).c[5],(w).c[7]); \
      _complex_mul((u).c[5],(v).c[3],(w).c[2]);\
      _complex_mul_assign((u).c[5],(v).c[4],(w).c[5]); \
      _complex_mul_assign((u).c[5],(v).c[5],(w).c[8]); \
      _complex_mul((u).c[6],(v).c[6],(w).c[0]);\
      _complex_mul_assign((u).c[6],(v).c[7],(w).c[3]); \
      _complex_mul_assign((u).c[6],(v).c[8],(w).c[6]); \
      _complex_mul((u).c[7],(v).c[6],(w).c[1]);\
      _complex_mul_assign((u).c[7],(v).c[7],(w).c[4]); \
      _complex_mul_assign((u).c[7],(v).c[8],(w).c[7]); \
      _complex_mul((u).c[8],(v).c[6],(w).c[2]);\
      _complex_mul_assign((u).c[8],(v).c[7],(w).c[5]); \
      _complex_mul_assign((u).c[8],(v).c[8],(w).c[8])

/* u=v*w^+ */
#define _suNg_times_suNg_dagger(u,v,w) \
      _complex_mul_star((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_star_assign((u).c[0],(v).c[1],(w).c[1]); \
      _complex_mul_star_assign((u).c[0],(v).c[2],(w).c[2]); \
      _complex_mul_star((u).c[1],(v).c[0],(w).c[3]);\
      _complex_mul_star_assign((u).c[1],(v).c[1],(w).c[4]); \
      _complex_mul_star_assign((u).c[1],(v).c[2],(w).c[5]); \
      _complex_mul_star((u).c[2],(v).c[0],(w).c[6]);\
      _complex_mul_star_assign((u).c[2],(v).c[1],(w).c[7]); \
      _complex_mul_star_assign((u).c[2],(v).c[2],(w).c[8]); \
      _complex_mul_star((u).c[3],(v).c[3],(w).c[0]);\
      _complex_mul_star_assign((u).c[3],(v).c[4],(w).c[1]); \
      _complex_mul_star_assign((u).c[3],(v).c[5],(w).c[2]); \
      _complex_mul_star((u).c[4],(v).c[3],(w).c[3]);\
      _complex_mul_star_assign((u).c[4],(v).c[4],(w).c[4]); \
      _complex_mul_star_assign((u).c[4],(v).c[5],(w).c[5]); \
      _complex_mul_star((u).c[5],(v).c[3],(w).c[6]);\
      _complex_mul_star_assign((u).c[5],(v).c[4],(w).c[7]); \
      _complex_mul_star_assign((u).c[5],(v).c[5],(w).c[8]); \
      _complex_mul_star((u).c[6],(v).c[6],(w).c[0]);\
      _complex_mul_star_assign((u).c[6],(v).c[7],(w).c[1]); \
      _complex_mul_star_assign((u).c[6],(v).c[8],(w).c[2]); \
      _complex_mul_star((u).c[7],(v).c[6],(w).c[3]);\
      _complex_mul_star_assign((u).c[7],(v).c[7],(w).c[4]); \
      _complex_mul_star_assign((u).c[7],(v).c[8],(w).c[5]); \
      _complex_mul_star((u).c[8],(v).c[6],(w).c[6]);\
      _complex_mul_star_assign((u).c[8],(v).c[7],(w).c[7]); \
      _complex_mul_star_assign((u).c[8],(v).c[8],(w).c[8])

/* u=v^+*w */
#define _suNg_dagger_times_suNg(u,v,w) \
      _complex_mul_star((u).c[0],(w).c[0],(v).c[0]);\
      _complex_mul_star_assign((u).c[0],(w).c[3],(v).c[3]); \
      _complex_mul_star_assign((u).c[0],(w).c[6],(v).c[6]); \
      _complex_mul_star((u).c[1],(w).c[1],(v).c[0]);\
      _complex_mul_star_assign((u).c[1],(w).c[4],(v).c[3]); \
      _complex_mul_star_assign((u).c[1],(w).c[7],(v).c[6]); \
      _complex_mul_star((u).c[2],(w).c[2],(v).c[0]);\
      _complex_mul_star_assign((u).c[2],(w).c[5],(v).c[3]); \
      _complex_mul_star_assign((u).c[2],(w).c[8],(v).c[6]); \
      _complex_mul_star((u).c[3],(w).c[0],(v).c[1]);\
      _complex_mul_star_assign((u).c[3],(w).c[3],(v).c[4]); \
      _complex_mul_star_assign((u).c[3],(w).c[6],(v).c[7]); \
      _complex_mul_star((u).c[4],(w).c[1],(v).c[1]);\
      _complex_mul_star_assign((u).c[4],(w).c[4],(v).c[4]); \
      _complex_mul_star_assign((u).c[4],(w).c[7],(v).c[7]); \
      _complex_mul_star((u).c[5],(w).c[2],(v).c[1]);\
      _complex_mul_star_assign((u).c[5],(w).c[5],(v).c[4]); \
      _complex_mul_star_assign((u).c[5],(w).c[8],(v).c[7]); \
      _complex_mul_star((u).c[6],(w).c[0],(v).c[2]);\
      _complex_mul_star_assign((u).c[6],(w).c[3],(v).c[5]); \
      _complex_mul_star_assign((u).c[6],(w).c[6],(v).c[8]); \
      _complex_mul_star((u).c[7],(w).c[1],(v).c[2]);\
      _complex_mul_star_assign((u).c[7],(w).c[4],(v).c[5]); \
      _complex_mul_star_assign((u).c[7],(w).c[7],(v).c[8]); \
      _complex_mul_star((u).c[8],(w).c[2],(v).c[2]);\
      _complex_mul_star_assign((u).c[8],(w).c[5],(v).c[5]); \
      _complex_mul_star_assign((u).c[8],(w).c[8],(v).c[8])

/* u=0 */
#define _suNg_zero(u) \
    _complex_0((u).c[0]);\
    _complex_0((u).c[1]);\
    _complex_0((u).c[2]);\
    _complex_0((u).c[3]);\
    _complex_0((u).c[4]);\
    _complex_0((u).c[5]);\
    _complex_0((u).c[6]);\
    _complex_0((u).c[7]);\
    _complex_0((u).c[8])

/* u=1 */
#define _suNg_unit(u) \
   _complex_1((u).c[0]);\
   _complex_0((u).c[1]);\
   _complex_0((u).c[2]);\
   _complex_0((u).c[3]);\
   _complex_1((u).c[4]);\
   _complex_0((u).c[5]);\
   _complex_0((u).c[6]);\
   _complex_0((u).c[7]);\
   _complex_1((u).c[8])

/* u=-v */
#define _suNg_minus(u,v) \
   _complex_minus((u).c[0],(v).c[0]);\
   _complex_minus((u).c[1],(v).c[1]);\
   _complex_minus((u).c[2],(v).c[2]);\
   _complex_minus((u).c[3],(v).c[3]);\
   _complex_minus((u).c[4],(v).c[4]);\
   _complex_minus((u).c[5],(v).c[5]);\
   _complex_minus((u).c[6],(v).c[6]);\
   _complex_minus((u).c[7],(v).c[7]);\
   _complex_minus((u).c[8],(v).c[8])

/* u=v */
#define _suNg_copy(u,v) \
   (u)=(v)

/* u=r*v (u,v mat, r real) */
#define _suNg_mul(u,r,v) \
   _complex_mulr((u).c[0],(r),(v).c[0]);\
   _complex_mulr((u).c[1],(r),(v).c[1]);\
   _complex_mulr((u).c[2],(r),(v).c[2]);\
   _complex_mulr((u).c[3],(r),(v).c[3]);\
   _complex_mulr((u).c[4],(r),(v).c[4]);\
   _complex_mulr((u).c[5],(r),(v).c[5]);\
   _complex_mulr((u).c[6],(r),(v).c[6]);\
   _complex_mulr((u).c[7],(r),(v).c[7]);\
   _complex_mulr((u).c[8],(r),(v).c[8])

/* u=r*v (u,v mat, r complex) */
#define _suNg_mulc(u,r,v) \
   _complex_mul((u).c[0],(r),(v).c[0]);\
   _complex_mul((u).c[1],(r),(v).c[1]);\
   _complex_mul((u).c[2],(r),(v).c[2]);\
   _complex_mul((u).c[3],(r),(v).c[3]);\
   _complex_mul((u).c[4],(r),(v).c[4]);\
   _complex_mul((u).c[5],(r),(v).c[5]);\
   _complex_mul((u).c[6],(r),(v).c[6]);\
   _complex_mul((u).c[7],(r),(v).c[7]);\
   _complex_mul((u).c[8],(r),(v).c[8])

/* u+=v */
#define _suNg_add_assign(u,v) \
   _complex_add_assign((u).c[0],(v).c[0]);\
   _complex_add_assign((u).c[1],(v).c[1]);\
   _complex_add_assign((u).c[2],(v).c[2]);\
   _complex_add_assign((u).c[3],(v).c[3]);\
   _complex_add_assign((u).c[4],(v).c[4]);\
   _complex_add_assign((u).c[5],(v).c[5]);\
   _complex_add_assign((u).c[6],(v).c[6]);\
   _complex_add_assign((u).c[7],(v).c[7]);\
   _complex_add_assign((u).c[8],(v).c[8])

/* u-=v */
#define _suNg_sub_assign(u,v) \
   _complex_sub_assign((u).c[0],(v).c[0]);\
   _complex_sub_assign((u).c[1],(v).c[1]);\
   _complex_sub_assign((u).c[2],(v).c[2]);\
   _complex_sub_assign((u).c[3],(v).c[3]);\
   _complex_sub_assign((u).c[4],(v).c[4]);\
   _complex_sub_assign((u).c[5],(v).c[5]);\
   _complex_sub_assign((u).c[6],(v).c[6]);\
   _complex_sub_assign((u).c[7],(v).c[7]);\
   _complex_sub_assign((u).c[8],(v).c[8])

/* k=| u |2 ) */
#define _suNg_sqnorm(k,u) \
   (k)=0.;\
   (k)+=_complex_prod_re((u).c[0],(u).c[0]); \
   (k)+=_complex_prod_re((u).c[1],(u).c[1]); \
   (k)+=_complex_prod_re((u).c[2],(u).c[2]); \
   (k)+=_complex_prod_re((u).c[3],(u).c[3]); \
   (k)+=_complex_prod_re((u).c[4],(u).c[4]); \
   (k)+=_complex_prod_re((u).c[5],(u).c[5]); \
   (k)+=_complex_prod_re((u).c[6],(u).c[6]); \
   (k)+=_complex_prod_re((u).c[7],(u).c[7]); \
   (k)+=_complex_prod_re((u).c[8],(u).c[8])

/* k=| 1 - u |2 ) */
#define _suNg_sqnorm_m1(k,u) \
   (k)=\
    +_complex_prod_m1_re((u).c[0],(u).c[0])\
    +_complex_prod_re((u).c[1],(u).c[1])\
    +_complex_prod_re((u).c[2],(u).c[2])\
    +_complex_prod_re((u).c[3],(u).c[3])\
    +_complex_prod_m1_re((u).c[4],(u).c[4])\
    +_complex_prod_re((u).c[5],(u).c[5])\
    +_complex_prod_re((u).c[6],(u).c[6])\
    +_complex_prod_re((u).c[7],(u).c[7])\
    +_complex_prod_m1_re((u).c[8],(u).c[8])

/* k=Re Tr (u) */
#define _suNg_trace_re(k,u) \
   (k)=_complex_re((u).c[0])+ \
       _complex_re((u).c[4])+ \
       _complex_re((u).c[8])

/* k=Im Tr (u) */
#define _suNg_trace_im(k,u) \
   (k)=_complex_im((u).c[0])+ \
       _complex_im((u).c[4])+ \
       _complex_im((u).c[8])

/* This fuction computes the hmc force matrix */
/* this fuction accumulates on the original matrix u */
#define _suNg_FMAT(u,s) \
   _complex_mul_star_assign((u).c[0],(s).c[0].c[0],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[0],(s).c[1].c[0],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[1],(s).c[0].c[0],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[1],(s).c[1].c[0],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[2],(s).c[0].c[0],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[2],(s).c[1].c[0],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[3],(s).c[0].c[1],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[3],(s).c[1].c[1],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[4],(s).c[0].c[1],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[4],(s).c[1].c[1],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[5],(s).c[0].c[1],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[5],(s).c[1].c[1],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[6],(s).c[0].c[2],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[6],(s).c[1].c[2],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[7],(s).c[0].c[2],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[7],(s).c[1].c[2],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[8],(s).c[0].c[2],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[8],(s).c[1].c[2],(s).c[3].c[2])

/* u=0 */
#define _suNg_FMAT_zero(u) \
    _complex_0((u).c[0]);\
    _complex_0((u).c[1]);\
    _complex_0((u).c[2]);\
    _complex_0((u).c[3]);\
    _complex_0((u).c[4]);\
    _complex_0((u).c[5]);\
    _complex_0((u).c[6]);\
    _complex_0((u).c[7]);\
    _complex_0((u).c[8])

/*******************************************************************************
*
* Macros for spinors
*
* Arguments are variables of type spinors
*
*******************************************************************************/

/*  r=0  (r spinor) */
#define _spinor_zero_g(r) \
  _vector_zero_g((r).c[0]); \
  _vector_zero_g((r).c[1]); \
  _vector_zero_g((r).c[2]); \
  _vector_zero_g((r).c[3])

/*  s=g5*r (r,s spinors, g5 matrix) */
#define _spinor_g5_g(s,r) \
  (s).c[0]=(r).c[0]; \
  (s).c[1]=(r).c[1]; \
  _vector_minus_g((s).c[2],(r).c[2]); \
  _vector_minus_g((s).c[3],(r).c[3])

/*  r=g5*r (r,s spinors, g5 matrix) */
#define _spinor_g5_assign_g(r) \
  _vector_minus_g((r).c[2],(r).c[2]); \
  _vector_minus_g((r).c[3],(r).c[3])

/*  s=-r (r,s spinors) */
#define _spinor_minus_g(s,r) \
  _vector_minus_g((s).c[0],(r).c[0]); \
  _vector_minus_g((s).c[1],(r).c[1]); \
  _vector_minus_g((s).c[2],(r).c[2]); \
  _vector_minus_g((s).c[3],(r).c[3])

/*  r=k*s (k real; r,s spinors) */
#define _spinor_mul_g(r,k,s) \
  _vector_mul_g((r).c[0],k,(s).c[0]); \
  _vector_mul_g((r).c[1],k,(s).c[1]); \
  _vector_mul_g((r).c[2],k,(s).c[2]); \
  _vector_mul_g((r).c[3],k,(s).c[3])

/*  r=z*s (z complex; r,s spinors) */
#define _spinor_mulc_g(r,z,s) \
  _vector_mulc_g((r).c[0],z,(s).c[0]); \
  _vector_mulc_g((r).c[1],z,(s).c[1]); \
  _vector_mulc_g((r).c[2],z,(s).c[2]); \
  _vector_mulc_g((r).c[3],z,(s).c[3])

/*  r+=z*s (z complex; r,s spinors) */
#define _spinor_mulc_add_assign_g(r,z,s) \
  _vector_mulc_add_assign_g((r).c[0],(z),(s).c[0]); \
  _vector_mulc_add_assign_g((r).c[1],(z),(s).c[1]); \
  _vector_mulc_add_assign_g((r).c[2],(z),(s).c[2]); \
  _vector_mulc_add_assign_g((r).c[3],(z),(s).c[3])

/*  r+=k*s (k real; r,s spinors) */
#define _spinor_mul_add_assign_g(r,k,s) \
  _vector_mul_add_assign_g((r).c[0],(k),(s).c[0]); \
  _vector_mul_add_assign_g((r).c[1],(k),(s).c[1]); \
  _vector_mul_add_assign_g((r).c[2],(k),(s).c[2]); \
  _vector_mul_add_assign_g((r).c[3],(k),(s).c[3])

/*  r=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */
#define _spinor_lc_g(r,k1,s1,k2,s2) \
  _vector_lc_g((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
  _vector_lc_g((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
  _vector_lc_g((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
  _vector_lc_g((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/*  r+=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */
#define _spinor_lc_add_assign_g(r,k1,s1,k2,s2) \
  _vector_lc_add_assign_g((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
  _vector_lc_add_assign_g((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
  _vector_lc_add_assign_g((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
  _vector_lc_add_assign_g((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/*  r=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */
#define _spinor_clc_g(r,z1,s1,z2,s2) \
  _vector_clc_g((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
  _vector_clc_g((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
  _vector_clc_g((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
  _vector_clc_g((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/*  r+=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */
#define _spinor_clc_add_assign_g(r,z1,s1,z2,s2) \
  _vector_clc_add_assign_g((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
  _vector_clc_add_assign_g((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
  _vector_clc_add_assign_g((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
  _vector_clc_add_assign_g((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/*  r=s1+s2 (r,s1,s2 spinors) */
#define _spinor_add_g(r,s1,s2) \
  _vector_add_g((r).c[0],(s1).c[0],(s2).c[0]); \
  _vector_add_g((r).c[1],(s1).c[1],(s2).c[1]); \
  _vector_add_g((r).c[2],(s1).c[2],(s2).c[2]); \
  _vector_add_g((r).c[3],(s1).c[3],(s2).c[3])

/*  r=s1-s2 (r,s1,s2 spinors) */
#define _spinor_sub_g(r,s1,s2) \
  _vector_sub_g((r).c[0],(s1).c[0],(s2).c[0]); \
  _vector_sub_g((r).c[1],(s1).c[1],(s2).c[1]); \
  _vector_sub_g((r).c[2],(s1).c[2],(s2).c[2]); \
  _vector_sub_g((r).c[3],(s1).c[3],(s2).c[3])

/*  r+=s (r,s spinors) */
#define _spinor_add_assign_g(r,s) \
  _vector_add_assign_g((r).c[0],(s).c[0]); \
  _vector_add_assign_g((r).c[1],(s).c[1]); \
  _vector_add_assign_g((r).c[2],(s).c[2]); \
  _vector_add_assign_g((r).c[3],(s).c[3])

/*  r-=s (r,s spinors) */
#define _spinor_sub_assign_g(r,s) \
  _vector_sub_assign_g((r).c[0],(s).c[0]); \
  _vector_sub_assign_g((r).c[1],(s).c[1]); \
  _vector_sub_assign_g((r).c[2],(s).c[2]); \
  _vector_sub_assign_g((r).c[3],(s).c[3])

/*  r+=i*s (r,s spinors) */
#define _spinor_i_add_assign_g(r,s) \
  _vector_i_add_assign_g((r).c[0],(s).c[0]); \
  _vector_i_add_assign_g((r).c[1],(s).c[1]); \
  _vector_i_add_assign_g((r).c[2],(s).c[2]); \
  _vector_i_add_assign_g((r).c[3],(s).c[3])

/*  r-=i*s (r,s spinors) */
#define _spinor_i_sub_assign_g(r,s) \
  _vector_i_sub_assign_g((r).c[0],(s).c[0]); \
  _vector_i_sub_assign_g((r).c[1],(s).c[1]); \
  _vector_i_sub_assign_g((r).c[2],(s).c[2]); \
  _vector_i_sub_assign_g((r).c[3],(s).c[3])

/* k=Real part of the scalar product r*s (r,s spinors) */
#define _spinor_prod_re_g(k,r,s) \
   do { \
      _vector_prod_re_g((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_re_g((k),(r).c[1],(s).c[1]); \
      _vector_prod_add_assign_re_g((k),(r).c[2],(s).c[2]); \
      _vector_prod_add_assign_re_g((k),(r).c[3],(s).c[3]); \
   } while(0) 

/* k=Im part of the scalar product r*s (r,s spinors) */
#define _spinor_prod_im_g(k,r,s) \
   do { \
      _vector_prod_im_g((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_im_g((k),(r).c[1],(s).c[1]); \
      _vector_prod_add_assign_im_g((k),(r).c[2],(s).c[2]); \
      _vector_prod_add_assign_im_g((k),(r).c[3],(s).c[3]); \
   } while(0) 

/* z=r*s (r,s spinors, z complex) */
#define _spinor_prod_g(z,r,s) \
   do { \
      (z).re=0.;(z).im=0.; \
      _vector_prod_assign_g((z),(r).c[0],(s).c[0]); \
      _vector_prod_assign_g((z),(r).c[1],(s).c[1]); \
      _vector_prod_assign_g((z),(r).c[2],(s).c[2]); \
      _vector_prod_assign_g((z),(r).c[3],(s).c[3]); \
   } while(0) 

/* z+=r*s (r,s spinors, z complex) */
#define _spinor_prod_assign_g(z,r,s) \
  _vector_prod_assign_g((z),(r).c[0],(s).c[0]); \
  _vector_prod_assign_g((z),(r).c[1],(s).c[1]); \
  _vector_prod_assign_g((z),(r).c[2],(s).c[2]); \
  _vector_prod_assign_g((z),(r).c[3],(s).c[3])

/* k=Real part of the scalar product (g5*r)*s (r,s spinors) */
#define _spinor_g5_prod_re_g(k,r,s) \
   do { \
      _vector_prod_re_g((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_re_g((k),(r).c[1],(s).c[1]);\
      _vector_prod_sub_assign_re_g((k),(r).c[2],(s).c[2]);\
      _vector_prod_sub_assign_re_g((k),(r).c[3],(s).c[3]);\
   } while(0) 

/* k=Imaginary part of the scalar product (g5*r)*s (r,s spinors) */
#define _spinor_g5_prod_im_g(k,r,s) \
   do { \
      _vector_prod_im_g((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_im_g((k),(r).c[1],(s).c[1]);\
      _vector_prod_sub_assign_im_g((k),(r).c[2],(s).c[2]);\
      _vector_prod_sub_assign_im_g((k),(r).c[3],(s).c[3]);\
   } while(0) 

/* r-=z*s (z complex; r,s spinors) */
#define _spinor_project_g(r,z,s) \
  _vector_project_g((r).c[0],z,(s).c[0]); \
  _vector_project_g((r).c[1],z,(s).c[1]); \
  _vector_project_g((r).c[2],z,(s).c[2]); \
  _vector_project_g((r).c[3],z,(s).c[3])

/* r=(1-g0)/2 * s (r,s spinors) */
#define _spinor_pminus_g(r,s) \
  _vector_add_g((r).c[0],(s).c[0],(s).c[2]); \
  _vector_add_g((r).c[1],(s).c[1],(s).c[3]); \
  _vector_mul_g((r).c[0],0.5,(r).c[0]); \
  _vector_mul_g((r).c[1],0.5,(r).c[1]); \
  (r).c[2] = (r).c[0]; \
  (r).c[3] = (r).c[1]

/* r=(1+g0)/2 * s (r,s spinors) */
#define _spinor_pplus_g(r,s) \
  _vector_sub_g((r).c[0],(s).c[0],(s).c[2]); \
  _vector_sub_g((r).c[1],(s).c[1],(s).c[3]); \
  _vector_mul_g((r).c[0],0.5,(r).c[0]); \
  _vector_mul_g((r).c[1],0.5,(r).c[1]); \
  _vector_mul_g((r).c[2],-1.,(r).c[0]); \
  _vector_mul_g((r).c[3],-1.,(r).c[1])

/* Read spinor field component from GPU memory */
/* (output) v = suNg_vector ; (input) in = suNg_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNg_read_spinor_flt_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*3)*(stride); \
      (v).c[0]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex_flt*)(in))[__iz]; \
   } while (0) 

#define _suNg_read_spinor_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*3)*(stride); \
      (v).c[0]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex*)(in))[__iz]; \
   } while (0) 

/* Write spinor field component to GPU memory */
/* (input) v = suNg_vector ; (output) out = suNg_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNg_write_spinor_flt_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*3)*(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[2]; \
   } while (0) 

#define _suNg_write_spinor_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*3)*(stride); \
      ((complex*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[2]; \
   } while (0) 

/* Read an suN algebra vector from GPU memory */
/* (output) v = suN_algebra_vector ; (input) in = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_av_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      (v).c[0]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((float*)(in))[__iz]; \
   } while (0) 

#define _suNg_av_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      (v).c[0]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN algebra vector to GPU memory */
/* (input) v = suN_algebra_vector ; (output) out = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_av_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      ((float*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7]; \
   } while (0) 

#define _suNg_av_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      ((double*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7]; \
   } while (0) 

/* Mul_add_assign on a suN algebra vector on GPU  */
/* (in/out) v = suN_algebra_vector* ; (input) in = suN_algebra_vector */
/* (input) iy = site ; (input) x = 0..3 direction; (input) r = real */
#define _algebra_vector_mul_add_assign_gpu_g_flt(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      ((float*)(v))[__iz]+=(in).c[0]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[1]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[2]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[3]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[4]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[5]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[6]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[7]*(r); \
   } while (0) 

#define _algebra_vector_mul_add_assign_gpu_g(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*8)*(stride); \
      ((double*)(v))[__iz]+=(in).c[0]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[1]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[2]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[3]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[4]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[5]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[6]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[7]*(r); \
   } while (0) 

/* Read an suN matrix from GPU memory */
/* (output) v = suN ; (input) in = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*18)*(stride); \
      (v).c[0].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[0].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].im=((float*)(in))[__iz]; \
   } while (0) 

#define _suNg_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*18)*(stride); \
      (v).c[0].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[0].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].im=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN matrix to GPU memory */
/* (input) v = suN ; (output) out = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*18)*(stride); \
      ((float*)(out))[__iz]=(v).c[0].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[0].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[8].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[8].im; \
   } while (0) 

#define _suNg_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*18)*(stride); \
      ((double*)(out))[__iz]=(v).c[0].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[0].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[8].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[8].im; \
   } while (0) 

/*******************************************************************************
*
* The following macros are the same for single and double precision types
*
* Depending on the macro, arguments are variables of type suN_vector and suN
* (or suN_vector_flt and suN_flt)
*
*******************************************************************************/

/* r=1 */
#define _vector_one_f(r) \
   _complex_1((r).c[0]); \
   _complex_1((r).c[1]); \
   _complex_1((r).c[2]); \
   _complex_1((r).c[3]); \
   _complex_1((r).c[4]); \
   _complex_1((r).c[5])

/* r=0 */
#define _vector_zero_f(r) \
   _complex_0((r).c[0]); \
   _complex_0((r).c[1]); \
   _complex_0((r).c[2]); \
   _complex_0((r).c[3]); \
   _complex_0((r).c[4]); \
   _complex_0((r).c[5])

/* r=-s */
#define _vector_minus_f(r,s) \
   _complex_minus((r).c[0],(s).c[0]); \
   _complex_minus((r).c[1],(s).c[1]); \
   _complex_minus((r).c[2],(s).c[2]); \
   _complex_minus((r).c[3],(s).c[3]); \
   _complex_minus((r).c[4],(s).c[4]); \
   _complex_minus((r).c[5],(s).c[5])

/* r= i*s */
#define _vector_i_plus_f(r,s) \
   _complex_i_plus((r).c[0],(s).c[0]); \
   _complex_i_plus((r).c[1],(s).c[1]); \
   _complex_i_plus((r).c[2],(s).c[2]); \
   _complex_i_plus((r).c[3],(s).c[3]); \
   _complex_i_plus((r).c[4],(s).c[4]); \
   _complex_i_plus((r).c[5],(s).c[5])

/* r=-i*s */
#define _vector_i_minus_f(r,s) \
   _complex_i_minus((r).c[0],(s).c[0]); \
   _complex_i_minus((r).c[1],(s).c[1]); \
   _complex_i_minus((r).c[2],(s).c[2]); \
   _complex_i_minus((r).c[3],(s).c[3]); \
   _complex_i_minus((r).c[4],(s).c[4]); \
   _complex_i_minus((r).c[5],(s).c[5])

/* r=k*s (k real) */
#define _vector_mul_f(r,k,s) \
   _complex_mulr((r).c[0],(k),(s).c[0]); \
   _complex_mulr((r).c[1],(k),(s).c[1]); \
   _complex_mulr((r).c[2],(k),(s).c[2]); \
   _complex_mulr((r).c[3],(k),(s).c[3]); \
   _complex_mulr((r).c[4],(k),(s).c[4]); \
   _complex_mulr((r).c[5],(k),(s).c[5])

/* r=z*s (z complex) */
#define _vector_mulc_f(r,z,s) \
   _complex_mul((r).c[0],(z),(s).c[0]); \
   _complex_mul((r).c[1],(z),(s).c[1]); \
   _complex_mul((r).c[2],(z),(s).c[2]); \
   _complex_mul((r).c[3],(z),(s).c[3]); \
   _complex_mul((r).c[4],(z),(s).c[4]); \
   _complex_mul((r).c[5],(z),(s).c[5])

/* r=(z^+)*s (z complex) */
#define _vector_mulc_star_f(r,z,s) \
   _complex_mul_star((r).c[0],(s).c[0],(z)); \
   _complex_mul_star((r).c[1],(s).c[1],(z)); \
   _complex_mul_star((r).c[2],(s).c[2],(z)); \
   _complex_mul_star((r).c[3],(s).c[3],(z)); \
   _complex_mul_star((r).c[4],(s).c[4],(z)); \
   _complex_mul_star((r).c[5],(s).c[5],(z))

/* r=s1+s2 */
#define _vector_add_f(r,s1,s2) \
   _complex_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_add((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_add((r).c[3],(s1).c[3],(s2).c[3]); \
   _complex_add((r).c[4],(s1).c[4],(s2).c[4]); \
   _complex_add((r).c[5],(s1).c[5],(s2).c[5])

/* r=s1-s2 */
#define _vector_sub_f(r,s1,s2) \
   _complex_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_sub((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_sub((r).c[3],(s1).c[3],(s2).c[3]); \
   _complex_sub((r).c[4],(s1).c[4],(s2).c[4]); \
   _complex_sub((r).c[5],(s1).c[5],(s2).c[5])

/* r=s1+i*s2 */
#define _vector_i_add_f(r,s1,s2) \
   _complex_i_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_add((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_i_add((r).c[3],(s1).c[3],(s2).c[3]); \
   _complex_i_add((r).c[4],(s1).c[4],(s2).c[4]); \
   _complex_i_add((r).c[5],(s1).c[5],(s2).c[5])

/* r=s1-i*s2 */
#define _vector_i_sub_f(r,s1,s2) \
   _complex_i_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_sub((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_i_sub((r).c[3],(s1).c[3],(s2).c[3]); \
   _complex_i_sub((r).c[4],(s1).c[4],(s2).c[4]); \
   _complex_i_sub((r).c[5],(s1).c[5],(s2).c[5])

/* r+=s */
#define _vector_add_assign_f(r,s) \
   _complex_add_assign((r).c[0],(s).c[0]); \
   _complex_add_assign((r).c[1],(s).c[1]); \
   _complex_add_assign((r).c[2],(s).c[2]); \
   _complex_add_assign((r).c[3],(s).c[3]); \
   _complex_add_assign((r).c[4],(s).c[4]); \
   _complex_add_assign((r).c[5],(s).c[5])

/* r-=s */
#define _vector_sub_assign_f(r,s) \
   _complex_sub_assign((r).c[0],(s).c[0]); \
   _complex_sub_assign((r).c[1],(s).c[1]); \
   _complex_sub_assign((r).c[2],(s).c[2]); \
   _complex_sub_assign((r).c[3],(s).c[3]); \
   _complex_sub_assign((r).c[4],(s).c[4]); \
   _complex_sub_assign((r).c[5],(s).c[5])

/* r+=i*s */
#define _vector_i_add_assign_f(r,s) \
   _complex_i_add_assign((r).c[0],(s).c[0]); \
   _complex_i_add_assign((r).c[1],(s).c[1]); \
   _complex_i_add_assign((r).c[2],(s).c[2]); \
   _complex_i_add_assign((r).c[3],(s).c[3]); \
   _complex_i_add_assign((r).c[4],(s).c[4]); \
   _complex_i_add_assign((r).c[5],(s).c[5])

/* r-=i*s */
#define _vector_i_sub_assign_f(r,s) \
   _complex_i_sub_assign((r).c[0],(s).c[0]); \
   _complex_i_sub_assign((r).c[1],(s).c[1]); \
   _complex_i_sub_assign((r).c[2],(s).c[2]); \
   _complex_i_sub_assign((r).c[3],(s).c[3]); \
   _complex_i_sub_assign((r).c[4],(s).c[4]); \
   _complex_i_sub_assign((r).c[5],(s).c[5])

/* k=Re(r^*s) */
#define _vector_prod_re_f(k,r,s) \
   (k)=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_re((r).c[3],(s).c[3]); \
   (k)+=_complex_prod_re((r).c[4],(s).c[4]); \
   (k)+=_complex_prod_re((r).c[5],(s).c[5])

/* k=Im(r*s) */
#define _vector_prod_im_f(k,r,s) \
   (k)=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_im((r).c[3],(s).c[3]); \
   (k)+=_complex_prod_im((r).c[4],(s).c[4]); \
   (k)+=_complex_prod_im((r).c[5],(s).c[5])

/* r+=z*s (z complex) */
#define _vector_mulc_add_assign_f(r,z,s) \
   _complex_mul_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_assign((r).c[2],(z),(s).c[2]); \
   _complex_mul_assign((r).c[3],(z),(s).c[3]); \
   _complex_mul_assign((r).c[4],(z),(s).c[4]); \
   _complex_mul_assign((r).c[5],(z),(s).c[5])

/* r+=k*s (k real) */
#define _vector_mul_add_assign_f(r,k,s) \
   _complex_mulr_assign((r).c[0],(k),(s).c[0]); \
   _complex_mulr_assign((r).c[1],(k),(s).c[1]); \
   _complex_mulr_assign((r).c[2],(k),(s).c[2]); \
   _complex_mulr_assign((r).c[3],(k),(s).c[3]); \
   _complex_mulr_assign((r).c[4],(k),(s).c[4]); \
   _complex_mulr_assign((r).c[5],(k),(s).c[5])

/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_f(r,k1,s1,k2,s2) \
   _complex_rlc((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
   _complex_rlc((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3]); \
   _complex_rlc((r).c[4],(k1),(s1).c[4],(k2),(s2).c[4]); \
   _complex_rlc((r).c[5],(k1),(s1).c[5],(k2),(s2).c[5])

/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_add_assign_f(r,k1,s1,k2,s2) \
   _complex_rlc_assign((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc_assign((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc_assign((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
   _complex_rlc_assign((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3]); \
   _complex_rlc_assign((r).c[4],(k1),(s1).c[4],(k2),(s2).c[4]); \
   _complex_rlc_assign((r).c[5],(k1),(s1).c[5],(k2),(s2).c[5])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_f(r,z1,s1,z2,s2) \
   _complex_clc((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
   _complex_clc((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3]); \
   _complex_clc((r).c[4],(z1),(s1).c[4],(z2),(s2).c[4]); \
   _complex_clc((r).c[5],(z1),(s1).c[5],(z2),(s2).c[5])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_add_assign_f(r,z1,s1,z2,s2) \
   _complex_clc_assign((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc_assign((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc_assign((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
   _complex_clc_assign((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3]); \
   _complex_clc_assign((r).c[4],(z1),(s1).c[4],(z2),(s2).c[4]); \
   _complex_clc_assign((r).c[5],(z1),(s1).c[5],(z2),(s2).c[5])

/* z+=r^*s (c complex) */
#define _vector_prod_assign_f(z,r,s) \
   _complex_prod_assign((z),(r).c[0],(s).c[0]); \
   _complex_prod_assign((z),(r).c[1],(s).c[1]); \
   _complex_prod_assign((z),(r).c[2],(s).c[2]); \
   _complex_prod_assign((z),(r).c[3],(s).c[3]); \
   _complex_prod_assign((z),(r).c[4],(s).c[4]); \
   _complex_prod_assign((z),(r).c[5],(s).c[5])

/* k+=Re(r^*s) */
#define _vector_prod_add_assign_re_f(k,r,s) \
   (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_re((r).c[3],(s).c[3]); \
   (k)+=_complex_prod_re((r).c[4],(s).c[4]); \
   (k)+=_complex_prod_re((r).c[5],(s).c[5])

/* k+=Im(r*s) */
#define _vector_prod_add_assign_im_f(k,r,s) \
   (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_im((r).c[3],(s).c[3]); \
   (k)+=_complex_prod_im((r).c[4],(s).c[4]); \
   (k)+=_complex_prod_im((r).c[5],(s).c[5])

/* k-=Re(r^*s) */
#define _vector_prod_sub_assign_re_f(k,r,s) \
   (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)-=_complex_prod_re((r).c[3],(s).c[3]); \
   (k)-=_complex_prod_re((r).c[4],(s).c[4]); \
   (k)-=_complex_prod_re((r).c[5],(s).c[5])

/* k-=Im(r*s) */
#define _vector_prod_sub_assign_im_f(k,r,s) \
   (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)-=_complex_prod_im((r).c[3],(s).c[3]); \
   (k)-=_complex_prod_im((r).c[4],(s).c[4]); \
   (k)-=_complex_prod_im((r).c[5],(s).c[5])

/* r-=z*s (z complex) */
#define _vector_project_f(r,z,s) \
   _complex_mul_sub_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_sub_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_sub_assign((r).c[2],(z),(s).c[2]); \
   _complex_mul_sub_assign((r).c[3],(z),(s).c[3]); \
   _complex_mul_sub_assign((r).c[4],(z),(s).c[4]); \
   _complex_mul_sub_assign((r).c[5],(z),(s).c[5])

/* SU(N) matrix u times SU(N) vector s */
/* r=u*s */
#define _suNf_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<6; ++_i){\
         _complex_mul((r).c[_i],(u).c[_k],(s).c[0]); ++_k;\
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[1]); ++_k;\
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[2]); ++_k;\
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[3]); ++_k;\
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[4]); ++_k;\
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[5]); ++_k;\
      }\
   } while(0) 

/* SU(N) matrix u^dagger times SU(N) vector s */
/* r=(u^dagger)*s */
#define _suNf_inverse_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<6; ++_i){\
         _complex_mul_star((r).c[_i],(s).c[0],(u).c[_k]);\
         _k+=6; _complex_mul_star_assign((r).c[_i],(s).c[1],(u).c[_k]);\
         _k+=6; _complex_mul_star_assign((r).c[_i],(s).c[2],(u).c[_k]);\
         _k+=6; _complex_mul_star_assign((r).c[_i],(s).c[3],(u).c[_k]);\
         _k+=6; _complex_mul_star_assign((r).c[_i],(s).c[4],(u).c[_k]);\
         _k+=6; _complex_mul_star_assign((r).c[_i],(s).c[5],(u).c[_k]);\
         _k-=29;\
      }\
   } while(0) 

/* u=0 */
#define _suNf_zero(u) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/*******************************************************************************
*
* Macros for SU(N) matrices
*
* Arguments are variables of type suN
*
*******************************************************************************/

/* u=v^dagger */
#define _suNf_dagger(u,v) \
   do { \
      int _i,_n=0,_k=0;\
      for (_i=0; _i<6; ++_i){\
         _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k+=6; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k+=6; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k+=6; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k+=6; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k+=6; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k-=29;\
      }\
   } while(0) 

/* u=v*w */
#define _suNf_times_suNf(u,v,w) \
   do { \
      int _i,_y,_n=0,_k=0,_l=0;\
      for (_i=0; _i<6; ++_i){\
         for (_y=0; _y<6; ++_y){\
            _complex_mul((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; _l+=6; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; _l+=6; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; _l+=6; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; _l+=6; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; _l+=6; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_n; _k-=5; _l-=29;\
         } _k+=6; _l=0;\
      }\
   } while(0) 

/* u=v*w^+ */
#define _suNf_times_suNf_dagger(u,v,w) \
   do { \
      int _i,_y,_n=0,_k=0,_l=0;\
      for (_i=0; _i<6; ++_i){\
         for (_y=0; _y<6; ++_y){\
            _complex_mul_star((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_n; _k-=5; ++_l;\
         } _k+=6; _l=0;\
      }\
   } while(0) 

/* u=v^+*w */
#define _suNf_dagger_times_suNf(u,v,w) \
   do { \
      int _i,_y,_n=0,_k=0,_l=0;\
      for (_i=0; _i<6; ++_i){\
         for (_y=0; _y<6; ++_y){\
            _k=_y; _l=_i;\
            _complex_mul_star((u).c[_n],(w).c[_k],(v).c[_l]);\
            _k+=6; _l+=6; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            _k+=6; _l+=6; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            _k+=6; _l+=6; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            _k+=6; _l+=6; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            _k+=6; _l+=6; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            ++_n;\
         }\
      }\
   } while(0) 

/* u=0 */
#define _suNf_zero(u) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=1 */
#define _suNf_unit(u) \
   do { \
      _suNf_zero((u));\
      _complex_1((u).c[0]);\
      _complex_1((u).c[7]);\
      _complex_1((u).c[14]);\
      _complex_1((u).c[21]);\
      _complex_1((u).c[28]);\
      _complex_1((u).c[35]);\
   } while(0) 

/* u=-v */
#define _suNf_minus(u,v) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=v */
#define _suNf_copy(u,v) \
   (u)=(v)

/* u=r*v (u,v mat, r real) */
#define _suNf_mul(u,r,v) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=r*v (u,v mat, r complex) */
#define _suNf_mulc(u,r,v) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u+=v */
#define _suNf_add_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u-=v */
#define _suNf_sub_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* k=| u |2 ) */
#define _suNf_sqnorm(k,u) \
   do { \
      int _i;\
      (k)=0.;\
      for (_i=0; _i<36; ){\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
      }\
   } while(0) 

/* k=| 1 - u |2 ) */
#define _suNf_sqnorm_m1(k,u) \
   (k)=\
    +_complex_prod_m1_re((u).c[0],(u).c[0])\
    +_complex_prod_re((u).c[1],(u).c[1])\
    +_complex_prod_re((u).c[2],(u).c[2])\
    +_complex_prod_re((u).c[3],(u).c[3])\
    +_complex_prod_re((u).c[4],(u).c[4])\
    +_complex_prod_re((u).c[5],(u).c[5])\
    +_complex_prod_re((u).c[6],(u).c[6])\
    +_complex_prod_m1_re((u).c[7],(u).c[7])\
    +_complex_prod_re((u).c[8],(u).c[8])\
    +_complex_prod_re((u).c[9],(u).c[9])\
    +_complex_prod_re((u).c[10],(u).c[10])\
    +_complex_prod_re((u).c[11],(u).c[11])\
    +_complex_prod_re((u).c[12],(u).c[12])\
    +_complex_prod_re((u).c[13],(u).c[13])\
    +_complex_prod_m1_re((u).c[14],(u).c[14])\
    +_complex_prod_re((u).c[15],(u).c[15])\
    +_complex_prod_re((u).c[16],(u).c[16])\
    +_complex_prod_re((u).c[17],(u).c[17])\
    +_complex_prod_re((u).c[18],(u).c[18])\
    +_complex_prod_re((u).c[19],(u).c[19])\
    +_complex_prod_re((u).c[20],(u).c[20])\
    +_complex_prod_m1_re((u).c[21],(u).c[21])\
    +_complex_prod_re((u).c[22],(u).c[22])\
    +_complex_prod_re((u).c[23],(u).c[23])\
    +_complex_prod_re((u).c[24],(u).c[24])\
    +_complex_prod_re((u).c[25],(u).c[25])\
    +_complex_prod_re((u).c[26],(u).c[26])\
    +_complex_prod_re((u).c[27],(u).c[27])\
    +_complex_prod_m1_re((u).c[28],(u).c[28])\
    +_complex_prod_re((u).c[29],(u).c[29])\
    +_complex_prod_re((u).c[30],(u).c[30])\
    +_complex_prod_re((u).c[31],(u).c[31])\
    +_complex_prod_re((u).c[32],(u).c[32])\
    +_complex_prod_re((u).c[33],(u).c[33])\
    +_complex_prod_re((u).c[34],(u).c[34])\
    +_complex_prod_m1_re((u).c[35],(u).c[35])

/* k=Re Tr (u) */
#define _suNf_trace_re(k,u) \
   (k)=_complex_re((u).c[0])+ \
       _complex_re((u).c[7])+ \
       _complex_re((u).c[14])+ \
       _complex_re((u).c[21])+ \
       _complex_re((u).c[28])+ \
       _complex_re((u).c[35])

/* k=Im Tr (u) */
#define _suNf_trace_im(k,u) \
   (k)=_complex_im((u).c[0])+ \
       _complex_im((u).c[7])+ \
       _complex_im((u).c[14])+ \
       _complex_im((u).c[21])+ \
       _complex_im((u).c[28])+ \
       _complex_im((u).c[35])

/* This fuction computes the hmc force matrix */
/* this fuction accumulates on the original matrix u */
#define _suNf_FMAT(u,s) \
   do { \
      int _i,_j,_n=0;\
      for (_i=0; _i<6; ++_i){\
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[0]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[0]); \
         ++_n; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[1]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[1]); \
         ++_n; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[2]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[2]); \
         ++_n; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[3]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[3]); \
         ++_n; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[4]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[4]); \
         ++_n; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[5]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[5]); \
         ++_n; \
      }\
   } while(0) 

#define _suNffull_multiply(a,b,c) _suNf_multiply(a,b,c)

#define _suNffull_mul(a,b,c) _suNf_mul(a,b,c)

#define _suNffull_add_assign(a,b) _suNf_add_assign(a,b)

#define _suNffull_sub_assign(a,b) _suNf_sub_assign(a,b)

#define _suNffull_times_suNffull(a,b,c) _suNf_times_suNf(a,b,c)

#define _suNffull_inverse_multiply(a,b,c) _suNf_inverse_multiply(a,b,c)

#define _suNffull_zero(a) _suNf_zero(a)

#define _suNf_expand(a,b) _suNf_copy(a,b)

/* u=0 */
#define _suNf_FMAT_zero(u) \
   do { \
      int _i;for (_i=0; _i<36; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/*******************************************************************************
*
* Macros for spinors
*
* Arguments are variables of type spinors
*
*******************************************************************************/

/*  r=1  (r spinor) */
#define _spinor_one_f(r) \
  _vector_one_f((r).c[0]); \
  _vector_one_f((r).c[1]); \
  _vector_one_f((r).c[2]); \
  _vector_one_f((r).c[3])

/*  r=0  (r spinor) */
#define _spinor_zero_f(r) \
  _vector_zero_f((r).c[0]); \
  _vector_zero_f((r).c[1]); \
  _vector_zero_f((r).c[2]); \
  _vector_zero_f((r).c[3])

/*  s=g5*r (r,s spinors, g5 matrix) */
#define _spinor_g5_f(s,r) \
  (s).c[0]=(r).c[0]; \
  (s).c[1]=(r).c[1]; \
  _vector_minus_f((s).c[2],(r).c[2]); \
  _vector_minus_f((s).c[3],(r).c[3])

/*  r=g5*r (r,s spinors, g5 matrix) */
#define _spinor_g5_assign_f(r) \
  _vector_minus_f((r).c[2],(r).c[2]); \
  _vector_minus_f((r).c[3],(r).c[3])

/*  s=-r (r,s spinors) */
#define _spinor_minus_f(s,r) \
  _vector_minus_f((s).c[0],(r).c[0]); \
  _vector_minus_f((s).c[1],(r).c[1]); \
  _vector_minus_f((s).c[2],(r).c[2]); \
  _vector_minus_f((s).c[3],(r).c[3])

/*  r=k*s (k real; r,s spinors) */
#define _spinor_mul_f(r,k,s) \
  _vector_mul_f((r).c[0],k,(s).c[0]); \
  _vector_mul_f((r).c[1],k,(s).c[1]); \
  _vector_mul_f((r).c[2],k,(s).c[2]); \
  _vector_mul_f((r).c[3],k,(s).c[3])

/*  r=z*s (z complex; r,s spinors) */
#define _spinor_mulc_f(r,z,s) \
  _vector_mulc_f((r).c[0],z,(s).c[0]); \
  _vector_mulc_f((r).c[1],z,(s).c[1]); \
  _vector_mulc_f((r).c[2],z,(s).c[2]); \
  _vector_mulc_f((r).c[3],z,(s).c[3])

/*  r+=z*s (z complex; r,s spinors) */
#define _spinor_mulc_add_assign_f(r,z,s) \
  _vector_mulc_add_assign_f((r).c[0],(z),(s).c[0]); \
  _vector_mulc_add_assign_f((r).c[1],(z),(s).c[1]); \
  _vector_mulc_add_assign_f((r).c[2],(z),(s).c[2]); \
  _vector_mulc_add_assign_f((r).c[3],(z),(s).c[3])

/*  r+=k*s (k real; r,s spinors) */
#define _spinor_mul_add_assign_f(r,k,s) \
  _vector_mul_add_assign_f((r).c[0],(k),(s).c[0]); \
  _vector_mul_add_assign_f((r).c[1],(k),(s).c[1]); \
  _vector_mul_add_assign_f((r).c[2],(k),(s).c[2]); \
  _vector_mul_add_assign_f((r).c[3],(k),(s).c[3])

/*  r=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */
#define _spinor_lc_f(r,k1,s1,k2,s2) \
  _vector_lc_f((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
  _vector_lc_f((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
  _vector_lc_f((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
  _vector_lc_f((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/*  r+=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */
#define _spinor_lc_add_assign_f(r,k1,s1,k2,s2) \
  _vector_lc_add_assign_f((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
  _vector_lc_add_assign_f((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
  _vector_lc_add_assign_f((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
  _vector_lc_add_assign_f((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/*  r=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */
#define _spinor_clc_f(r,z1,s1,z2,s2) \
  _vector_clc_f((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
  _vector_clc_f((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
  _vector_clc_f((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
  _vector_clc_f((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/*  r+=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */
#define _spinor_clc_add_assign_f(r,z1,s1,z2,s2) \
  _vector_clc_add_assign_f((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
  _vector_clc_add_assign_f((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
  _vector_clc_add_assign_f((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
  _vector_clc_add_assign_f((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/*  r=s1+s2 (r,s1,s2 spinors) */
#define _spinor_add_f(r,s1,s2) \
  _vector_add_f((r).c[0],(s1).c[0],(s2).c[0]); \
  _vector_add_f((r).c[1],(s1).c[1],(s2).c[1]); \
  _vector_add_f((r).c[2],(s1).c[2],(s2).c[2]); \
  _vector_add_f((r).c[3],(s1).c[3],(s2).c[3])

/*  r=s1-s2 (r,s1,s2 spinors) */
#define _spinor_sub_f(r,s1,s2) \
  _vector_sub_f((r).c[0],(s1).c[0],(s2).c[0]); \
  _vector_sub_f((r).c[1],(s1).c[1],(s2).c[1]); \
  _vector_sub_f((r).c[2],(s1).c[2],(s2).c[2]); \
  _vector_sub_f((r).c[3],(s1).c[3],(s2).c[3])

/*  r+=s (r,s spinors) */
#define _spinor_add_assign_f(r,s) \
  _vector_add_assign_f((r).c[0],(s).c[0]); \
  _vector_add_assign_f((r).c[1],(s).c[1]); \
  _vector_add_assign_f((r).c[2],(s).c[2]); \
  _vector_add_assign_f((r).c[3],(s).c[3])

/*  r-=s (r,s spinors) */
#define _spinor_sub_assign_f(r,s) \
  _vector_sub_assign_f((r).c[0],(s).c[0]); \
  _vector_sub_assign_f((r).c[1],(s).c[1]); \
  _vector_sub_assign_f((r).c[2],(s).c[2]); \
  _vector_sub_assign_f((r).c[3],(s).c[3])

/*  r+=i*s (r,s spinors) */
#define _spinor_i_add_assign_f(r,s) \
  _vector_i_add_assign_f((r).c[0],(s).c[0]); \
  _vector_i_add_assign_f((r).c[1],(s).c[1]); \
  _vector_i_add_assign_f((r).c[2],(s).c[2]); \
  _vector_i_add_assign_f((r).c[3],(s).c[3])

/*  r-=i*s (r,s spinors) */
#define _spinor_i_sub_assign_f(r,s) \
  _vector_i_sub_assign_f((r).c[0],(s).c[0]); \
  _vector_i_sub_assign_f((r).c[1],(s).c[1]); \
  _vector_i_sub_assign_f((r).c[2],(s).c[2]); \
  _vector_i_sub_assign_f((r).c[3],(s).c[3])

/* k=Real part of the scalar product r*s (r,s spinors) */
#define _spinor_prod_re_f(k,r,s) \
   do { \
      _vector_prod_re_f((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_re_f((k),(r).c[1],(s).c[1]); \
      _vector_prod_add_assign_re_f((k),(r).c[2],(s).c[2]); \
      _vector_prod_add_assign_re_f((k),(r).c[3],(s).c[3]); \
   } while(0) 

/* k=Im part of the scalar product r*s (r,s spinors) */
#define _spinor_prod_im_f(k,r,s) \
   do { \
      _vector_prod_im_f((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_im_f((k),(r).c[1],(s).c[1]); \
      _vector_prod_add_assign_im_f((k),(r).c[2],(s).c[2]); \
      _vector_prod_add_assign_im_f((k),(r).c[3],(s).c[3]); \
   } while(0) 

/* z=r*s (r,s spinors, z complex) */
#define _spinor_prod_f(z,r,s) \
   do { \
      (z).re=0.;(z).im=0.; \
      _vector_prod_assign_f((z),(r).c[0],(s).c[0]); \
      _vector_prod_assign_f((z),(r).c[1],(s).c[1]); \
      _vector_prod_assign_f((z),(r).c[2],(s).c[2]); \
      _vector_prod_assign_f((z),(r).c[3],(s).c[3]); \
   } while(0) 

/* z+=r*s (r,s spinors, z complex) */
#define _spinor_prod_assign_f(z,r,s) \
  _vector_prod_assign_f((z),(r).c[0],(s).c[0]); \
  _vector_prod_assign_f((z),(r).c[1],(s).c[1]); \
  _vector_prod_assign_f((z),(r).c[2],(s).c[2]); \
  _vector_prod_assign_f((z),(r).c[3],(s).c[3])

/* k=Real part of the scalar product (g5*r)*s (r,s spinors) */
#define _spinor_g5_prod_re_f(k,r,s) \
   do { \
      _vector_prod_re_f((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_re_f((k),(r).c[1],(s).c[1]);\
      _vector_prod_sub_assign_re_f((k),(r).c[2],(s).c[2]);\
      _vector_prod_sub_assign_re_f((k),(r).c[3],(s).c[3]);\
   } while(0) 

/* k=Imaginary part of the scalar product (g5*r)*s (r,s spinors) */
#define _spinor_g5_prod_im_f(k,r,s) \
   do { \
      _vector_prod_im_f((k),(r).c[0],(s).c[0]);\
      _vector_prod_add_assign_im_f((k),(r).c[1],(s).c[1]);\
      _vector_prod_sub_assign_im_f((k),(r).c[2],(s).c[2]);\
      _vector_prod_sub_assign_im_f((k),(r).c[3],(s).c[3]);\
   } while(0) 

/* r-=z*s (z complex; r,s spinors) */
#define _spinor_project_f(r,z,s) \
  _vector_project_f((r).c[0],z,(s).c[0]); \
  _vector_project_f((r).c[1],z,(s).c[1]); \
  _vector_project_f((r).c[2],z,(s).c[2]); \
  _vector_project_f((r).c[3],z,(s).c[3])

/* r=(1-g0)/2 * s (r,s spinors) */
#define _spinor_pminus_f(r,s) \
  _vector_add_f((r).c[0],(s).c[0],(s).c[2]); \
  _vector_add_f((r).c[1],(s).c[1],(s).c[3]); \
  _vector_mul_f((r).c[0],0.5,(r).c[0]); \
  _vector_mul_f((r).c[1],0.5,(r).c[1]); \
  (r).c[2] = (r).c[0]; \
  (r).c[3] = (r).c[1]

/* r=(1+g0)/2 * s (r,s spinors) */
#define _spinor_pplus_f(r,s) \
  _vector_sub_f((r).c[0],(s).c[0],(s).c[2]); \
  _vector_sub_f((r).c[1],(s).c[1],(s).c[3]); \
  _vector_mul_f((r).c[0],0.5,(r).c[0]); \
  _vector_mul_f((r).c[1],0.5,(r).c[1]); \
  _vector_mul_f((r).c[2],-1.,(r).c[0]); \
  _vector_mul_f((r).c[3],-1.,(r).c[1])

/* Read spinor field component from GPU memory */
/* (output) v = suNf_vector ; (input) in = suNf_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNf_read_spinor_flt_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*6)*(stride); \
      (v).c[0]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((complex_flt*)(in))[__iz]; \
   } while (0) 

#define _suNf_read_spinor_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*6)*(stride); \
      (v).c[0]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((complex*)(in))[__iz]; \
   } while (0) 

/* Write spinor field component to GPU memory */
/* (input) v = suNf_vector ; (output) out = suNf_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNf_write_spinor_flt_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*6)*(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[5]; \
   } while (0) 

#define _suNf_write_spinor_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*6)*(stride); \
      ((complex*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[5]; \
   } while (0) 

/* Read an suN algebra vector from GPU memory */
/* (output) v = suN_algebra_vector ; (input) in = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_av_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      (v).c[0]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[8]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[9]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[10]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[11]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[12]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[13]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[14]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[15]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[16]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[17]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[18]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[19]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[20]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[21]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[22]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[23]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[24]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[25]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[26]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[27]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[28]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[29]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[30]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[31]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[32]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[33]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[34]=((float*)(in))[__iz]; \
   } while (0) 

#define _suNf_av_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      (v).c[0]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[8]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[9]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[10]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[11]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[12]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[13]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[14]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[15]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[16]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[17]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[18]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[19]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[20]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[21]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[22]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[23]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[24]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[25]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[26]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[27]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[28]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[29]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[30]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[31]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[32]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[33]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[34]=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN algebra vector to GPU memory */
/* (input) v = suN_algebra_vector ; (output) out = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_av_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      ((float*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[8]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[9]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[10]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[11]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[12]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[13]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[14]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[15]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[16]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[17]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[18]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[19]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[20]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[21]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[22]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[23]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[24]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[25]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[26]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[27]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[28]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[29]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[30]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[31]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[32]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[33]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[34]; \
   } while (0) 

#define _suNf_av_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      ((double*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[8]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[9]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[10]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[11]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[12]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[13]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[14]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[15]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[16]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[17]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[18]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[19]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[20]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[21]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[22]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[23]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[24]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[25]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[26]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[27]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[28]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[29]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[30]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[31]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[32]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[33]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[34]; \
   } while (0) 

/* Mul_add_assign on a suN algebra vector on GPU  */
/* (in/out) v = suN_algebra_vector* ; (input) in = suN_algebra_vector */
/* (input) iy = site ; (input) x = 0..3 direction; (input) r = real */
#define _algebra_vector_mul_add_assign_gpu_f_flt(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      ((float*)(v))[__iz]+=(in).c[0]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[1]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[2]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[3]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[4]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[5]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[6]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[7]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[8]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[9]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[10]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[11]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[12]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[13]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[14]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[15]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[16]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[17]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[18]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[19]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[20]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[21]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[22]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[23]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[24]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[25]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[26]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[27]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[28]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[29]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[30]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[31]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[32]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[33]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[34]*(r); \
   } while (0) 

#define _algebra_vector_mul_add_assign_gpu_f(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*35)*(stride); \
      ((double*)(v))[__iz]+=(in).c[0]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[1]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[2]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[3]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[4]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[5]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[6]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[7]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[8]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[9]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[10]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[11]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[12]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[13]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[14]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[15]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[16]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[17]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[18]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[19]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[20]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[21]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[22]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[23]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[24]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[25]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[26]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[27]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[28]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[29]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[30]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[31]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[32]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[33]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[34]*(r); \
   } while (0) 

/* Read an suN matrix from GPU memory */
/* (output) v = suN ; (input) in = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*72)*(stride); \
      (v).c[0].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[0].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[9].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[9].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[10].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[10].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[11].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[11].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[12].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[12].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[13].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[13].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[14].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[14].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[15].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[15].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[16].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[16].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[17].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[17].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[18].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[18].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[19].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[19].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[20].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[20].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[21].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[21].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[22].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[22].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[23].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[23].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[24].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[24].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[25].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[25].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[26].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[26].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[27].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[27].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[28].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[28].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[29].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[29].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[30].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[30].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[31].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[31].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[32].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[32].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[33].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[33].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[34].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[34].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[35].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[35].im=((float*)(in))[__iz]; \
   } while (0) 

#define _suNf_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*72)*(stride); \
      (v).c[0].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[0].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[1].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[2].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[3].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[4].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[5].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[6].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[7].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[8].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[9].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[9].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[10].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[10].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[11].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[11].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[12].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[12].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[13].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[13].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[14].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[14].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[15].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[15].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[16].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[16].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[17].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[17].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[18].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[18].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[19].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[19].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[20].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[20].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[21].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[21].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[22].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[22].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[23].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[23].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[24].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[24].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[25].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[25].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[26].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[26].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[27].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[27].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[28].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[28].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[29].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[29].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[30].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[30].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[31].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[31].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[32].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[32].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[33].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[33].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[34].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[34].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[35].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[35].im=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN matrix to GPU memory */
/* (input) v = suN ; (output) out = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*72)*(stride); \
      ((float*)(out))[__iz]=(v).c[0].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[0].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[1].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[2].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[3].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[4].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[5].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[6].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[7].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[8].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[8].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[9].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[9].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[10].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[10].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[11].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[11].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[12].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[12].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[13].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[13].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[14].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[14].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[15].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[15].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[16].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[16].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[17].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[17].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[18].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[18].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[19].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[19].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[20].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[20].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[21].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[21].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[22].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[22].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[23].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[23].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[24].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[24].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[25].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[25].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[26].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[26].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[27].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[27].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[28].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[28].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[29].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[29].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[30].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[30].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[31].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[31].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[32].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[32].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[33].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[33].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[34].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[34].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[35].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[35].im; \
   } while (0) 

#define _suNf_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*72)*(stride); \
      ((double*)(out))[__iz]=(v).c[0].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[0].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[1].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[2].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[3].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[4].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[5].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[6].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[7].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[8].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[8].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[9].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[9].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[10].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[10].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[11].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[11].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[12].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[12].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[13].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[13].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[14].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[14].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[15].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[15].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[16].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[16].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[17].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[17].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[18].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[18].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[19].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[19].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[20].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[20].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[21].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[21].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[22].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[22].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[23].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[23].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[24].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[24].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[25].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[25].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[26].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[26].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[27].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[27].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[28].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[28].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[29].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[29].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[30].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[30].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[31].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[31].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[32].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[32].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[33].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[33].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[34].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[34].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[35].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[35].im; \
   } while (0) 


#endif
