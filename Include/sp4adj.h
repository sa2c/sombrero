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
   _complex_0((r).c[2]); \
   _complex_0((r).c[3])

/* r=-s */
#define _vector_minus_g(r,s) \
   _complex_minus((r).c[0],(s).c[0]); \
   _complex_minus((r).c[1],(s).c[1]); \
   _complex_minus((r).c[2],(s).c[2]); \
   _complex_minus((r).c[3],(s).c[3])

/* r= i*s */
#define _vector_i_plus_g(r,s) \
   _complex_i_plus((r).c[0],(s).c[0]); \
   _complex_i_plus((r).c[1],(s).c[1]); \
   _complex_i_plus((r).c[2],(s).c[2]); \
   _complex_i_plus((r).c[3],(s).c[3])

/* r=-i*s */
#define _vector_i_minus_g(r,s) \
   _complex_i_minus((r).c[0],(s).c[0]); \
   _complex_i_minus((r).c[1],(s).c[1]); \
   _complex_i_minus((r).c[2],(s).c[2]); \
   _complex_i_minus((r).c[3],(s).c[3])

/* r=k*s (k real) */
#define _vector_mul_g(r,k,s) \
   _complex_mulr((r).c[0],(k),(s).c[0]); \
   _complex_mulr((r).c[1],(k),(s).c[1]); \
   _complex_mulr((r).c[2],(k),(s).c[2]); \
   _complex_mulr((r).c[3],(k),(s).c[3])

/* r=z*s (z complex) */
#define _vector_mulc_g(r,z,s) \
   _complex_mul((r).c[0],(z),(s).c[0]); \
   _complex_mul((r).c[1],(z),(s).c[1]); \
   _complex_mul((r).c[2],(z),(s).c[2]); \
   _complex_mul((r).c[3],(z),(s).c[3])

/* r=(z^+)*s (z complex) */
#define _vector_mulc_star_g(r,z,s) \
   _complex_mul_star((r).c[0],(s).c[0],(z)); \
   _complex_mul_star((r).c[1],(s).c[1],(z)); \
   _complex_mul_star((r).c[2],(s).c[2],(z)); \
   _complex_mul_star((r).c[3],(s).c[3],(z))

/* r=s1+s2 */
#define _vector_add_g(r,s1,s2) \
   _complex_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_add((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_add((r).c[3],(s1).c[3],(s2).c[3])

/* r=s1-s2 */
#define _vector_sub_g(r,s1,s2) \
   _complex_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_sub((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_sub((r).c[3],(s1).c[3],(s2).c[3])

/* r=s1+i*s2 */
#define _vector_i_add_g(r,s1,s2) \
   _complex_i_add((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_add((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_add((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_i_add((r).c[3],(s1).c[3],(s2).c[3])

/* r=s1-i*s2 */
#define _vector_i_sub_g(r,s1,s2) \
   _complex_i_sub((r).c[0],(s1).c[0],(s2).c[0]); \
   _complex_i_sub((r).c[1],(s1).c[1],(s2).c[1]); \
   _complex_i_sub((r).c[2],(s1).c[2],(s2).c[2]); \
   _complex_i_sub((r).c[3],(s1).c[3],(s2).c[3])

/* r+=s */
#define _vector_add_assign_g(r,s) \
   _complex_add_assign((r).c[0],(s).c[0]); \
   _complex_add_assign((r).c[1],(s).c[1]); \
   _complex_add_assign((r).c[2],(s).c[2]); \
   _complex_add_assign((r).c[3],(s).c[3])

/* r-=s */
#define _vector_sub_assign_g(r,s) \
   _complex_sub_assign((r).c[0],(s).c[0]); \
   _complex_sub_assign((r).c[1],(s).c[1]); \
   _complex_sub_assign((r).c[2],(s).c[2]); \
   _complex_sub_assign((r).c[3],(s).c[3])

/* r+=i*s */
#define _vector_i_add_assign_g(r,s) \
   _complex_i_add_assign((r).c[0],(s).c[0]); \
   _complex_i_add_assign((r).c[1],(s).c[1]); \
   _complex_i_add_assign((r).c[2],(s).c[2]); \
   _complex_i_add_assign((r).c[3],(s).c[3])

/* r-=i*s */
#define _vector_i_sub_assign_g(r,s) \
   _complex_i_sub_assign((r).c[0],(s).c[0]); \
   _complex_i_sub_assign((r).c[1],(s).c[1]); \
   _complex_i_sub_assign((r).c[2],(s).c[2]); \
   _complex_i_sub_assign((r).c[3],(s).c[3])

/* k=Re(r^*s) */
#define _vector_prod_re_g(k,r,s) \
   (k)=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_re((r).c[3],(s).c[3])

/* k=Im(r*s) */
#define _vector_prod_im_g(k,r,s) \
   (k)=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_im((r).c[3],(s).c[3])

/* r+=z*s (z complex) */
#define _vector_mulc_add_assign_g(r,z,s) \
   _complex_mul_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_assign((r).c[2],(z),(s).c[2]); \
   _complex_mul_assign((r).c[3],(z),(s).c[3])

/* r+=k*s (k real) */
#define _vector_mul_add_assign_g(r,k,s) \
   _complex_mulr_assign((r).c[0],(k),(s).c[0]); \
   _complex_mulr_assign((r).c[1],(k),(s).c[1]); \
   _complex_mulr_assign((r).c[2],(k),(s).c[2]); \
   _complex_mulr_assign((r).c[3],(k),(s).c[3])

/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_g(r,k1,s1,k2,s2) \
   _complex_rlc((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
   _complex_rlc((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_add_assign_g(r,k1,s1,k2,s2) \
   _complex_rlc_assign((r).c[0],(k1),(s1).c[0],(k2),(s2).c[0]); \
   _complex_rlc_assign((r).c[1],(k1),(s1).c[1],(k2),(s2).c[1]); \
   _complex_rlc_assign((r).c[2],(k1),(s1).c[2],(k2),(s2).c[2]); \
   _complex_rlc_assign((r).c[3],(k1),(s1).c[3],(k2),(s2).c[3])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_g(r,z1,s1,z2,s2) \
   _complex_clc((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
   _complex_clc((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_add_assign_g(r,z1,s1,z2,s2) \
   _complex_clc_assign((r).c[0],(z1),(s1).c[0],(z2),(s2).c[0]); \
   _complex_clc_assign((r).c[1],(z1),(s1).c[1],(z2),(s2).c[1]); \
   _complex_clc_assign((r).c[2],(z1),(s1).c[2],(z2),(s2).c[2]); \
   _complex_clc_assign((r).c[3],(z1),(s1).c[3],(z2),(s2).c[3])

/* z+=r^*s (c complex) */
#define _vector_prod_assign_g(z,r,s) \
   _complex_prod_assign((z),(r).c[0],(s).c[0]); \
   _complex_prod_assign((z),(r).c[1],(s).c[1]); \
   _complex_prod_assign((z),(r).c[2],(s).c[2]); \
   _complex_prod_assign((z),(r).c[3],(s).c[3])

/* k+=Re(r^*s) */
#define _vector_prod_add_assign_re_g(k,r,s) \
   (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_re((r).c[3],(s).c[3])

/* k+=Im(r*s) */
#define _vector_prod_add_assign_im_g(k,r,s) \
   (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)+=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)+=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)+=_complex_prod_im((r).c[3],(s).c[3])

/* k-=Re(r^*s) */
#define _vector_prod_sub_assign_re_g(k,r,s) \
   (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_re((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_re((r).c[2],(s).c[2]); \
   (k)-=_complex_prod_re((r).c[3],(s).c[3])

/* k-=Im(r*s) */
#define _vector_prod_sub_assign_im_g(k,r,s) \
   (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
   (k)-=_complex_prod_im((r).c[1],(s).c[1]); \
   (k)-=_complex_prod_im((r).c[2],(s).c[2]); \
   (k)-=_complex_prod_im((r).c[3],(s).c[3])

/* r-=z*s (z complex) */
#define _vector_project_g(r,z,s) \
   _complex_mul_sub_assign((r).c[0],(z),(s).c[0]); \
   _complex_mul_sub_assign((r).c[1],(z),(s).c[1]); \
   _complex_mul_sub_assign((r).c[2],(z),(s).c[2]); \
   _complex_mul_sub_assign((r).c[3],(z),(s).c[3])

/* SP(N) matrix u times SP(N) vector s */
/* r=u*s */
#define _suNg_multiply(r,u,s) \
      _complex_mul((r).c[0],(u).c[0],(s).c[0]);\
      _complex_mul_assign((r).c[0],(u).c[1],(s).c[1]); \
      _complex_mul_assign((r).c[0],(u).c[2],(s).c[2]); \
      _complex_mul_assign((r).c[0],(u).c[3],(s).c[3]); \
      _complex_mul((r).c[1],(u).c[4],(s).c[0]);\
      _complex_mul_assign((r).c[1],(u).c[5],(s).c[1]); \
      _complex_mul_assign((r).c[1],(u).c[6],(s).c[2]); \
      _complex_mul_assign((r).c[1],(u).c[7],(s).c[3]); \
      _complex_minus_mul_star((r).c[2],(s).c[0],(u).c[2]);\
      _complex_mul_star_massign((r).c[2],(s).c[1],(u).c[3]);\
      _complex_mul_star_passign((r).c[2],(s).c[2],(u).c[0]);\
      _complex_mul_star_passign((r).c[2],(s).c[3],(u).c[1]);\
      _complex_minus_mul_star((r).c[3],(s).c[0],(u).c[6]);\
      _complex_mul_star_massign((r).c[3],(s).c[1],(u).c[7]);\
      _complex_mul_star_passign((r).c[3],(s).c[2],(u).c[4]);\
      _complex_mul_star_passign((r).c[3],(s).c[3],(u).c[5]);\


/* SP(N) matrix u^dagger times SP(N) vector s */
/* r=(u^dagger)*s */
#define _suNg_inverse_multiply(r,u,s) \
      _complex_mul_star((r).c[0],(s).c[0],(u).c[0]);\
      _complex_mul_star_passign((r).c[0],(s).c[1],(u).c[4]); \
      _complex_mul_massign((r).c[0],(s).c[2],(u).c[2]); \
      _complex_mul_massign((r).c[0],(s).c[3],(u).c[6]); \
      _complex_mul_star((r).c[1],(s).c[0],(u).c[1]);\
      _complex_mul_star_passign((r).c[1],(s).c[1],(u).c[5]); \
      _complex_mul_massign((r).c[1],(s).c[2],(u).c[3]); \
      _complex_mul_massign((r).c[1],(s).c[3],(u).c[7]); \
      _complex_mul_star((r).c[2],(s).c[0],(u).c[2]);\
      _complex_mul_star_passign((r).c[2],(s).c[1],(u).c[6]); \
      _complex_mul_passign((r).c[2],(s).c[2],(u).c[0]); \
      _complex_mul_passign((r).c[2],(s).c[3],(u).c[4]); \
      _complex_mul_star((r).c[3],(s).c[0],(u).c[3]);\
      _complex_mul_star_passign((r).c[3],(s).c[1],(u).c[7]); \
      _complex_mul_passign((r).c[3],(s).c[2],(u).c[1]); \
      _complex_mul_passign((r).c[3],(s).c[3],(u).c[5]); \


/* SU(N) matrix u times SU(N) vector s */
/* r=u*s */
#define _suNgfull_multiply(r,u,s) \
      _complex_mul((r).c[0],(u).c[0],(s).c[0]);\
      _complex_mul_assign((r).c[0],(u).c[1],(s).c[1]); \
      _complex_mul_assign((r).c[0],(u).c[2],(s).c[2]); \
      _complex_mul_assign((r).c[0],(u).c[3],(s).c[3]); \
      _complex_mul((r).c[1],(u).c[4],(s).c[0]);\
      _complex_mul_assign((r).c[1],(u).c[5],(s).c[1]); \
      _complex_mul_assign((r).c[1],(u).c[6],(s).c[2]); \
      _complex_mul_assign((r).c[1],(u).c[7],(s).c[3]); \
      _complex_mul((r).c[2],(u).c[8],(s).c[0]);\
      _complex_mul_assign((r).c[2],(u).c[9],(s).c[1]); \
      _complex_mul_assign((r).c[2],(u).c[10],(s).c[2]); \
      _complex_mul_assign((r).c[2],(u).c[11],(s).c[3]); \
      _complex_mul((r).c[3],(u).c[12],(s).c[0]);\
      _complex_mul_assign((r).c[3],(u).c[13],(s).c[1]); \
      _complex_mul_assign((r).c[3],(u).c[14],(s).c[2]); \
      _complex_mul_assign((r).c[3],(u).c[15],(s).c[3])

/* SU(N) matrix u^dagger times SU(N) vector s */
/* r=(u^dagger)*s */
#define _suNgfull_inverse_multiply(r,u,s) \
      _complex_mul_star((r).c[0],(s).c[0],(u).c[0]);\
      _complex_mul_star_assign((r).c[0],(s).c[1],(u).c[4]); \
      _complex_mul_star_assign((r).c[0],(s).c[2],(u).c[8]); \
      _complex_mul_star_assign((r).c[0],(s).c[3],(u).c[12]); \
      _complex_mul_star((r).c[1],(s).c[0],(u).c[1]);\
      _complex_mul_star_assign((r).c[1],(s).c[1],(u).c[5]); \
      _complex_mul_star_assign((r).c[1],(s).c[2],(u).c[9]); \
      _complex_mul_star_assign((r).c[1],(s).c[3],(u).c[13]); \
      _complex_mul_star((r).c[2],(s).c[0],(u).c[2]);\
      _complex_mul_star_assign((r).c[2],(s).c[1],(u).c[6]); \
      _complex_mul_star_assign((r).c[2],(s).c[2],(u).c[10]); \
      _complex_mul_star_assign((r).c[2],(s).c[3],(u).c[14]); \
      _complex_mul_star((r).c[3],(s).c[0],(u).c[3]);\
      _complex_mul_star_assign((r).c[3],(s).c[1],(u).c[7]); \
      _complex_mul_star_assign((r).c[3],(s).c[2],(u).c[11]); \
      _complex_mul_star_assign((r).c[3],(s).c[3],(u).c[15])

/* r+=s */
#define _algebra_vector_add_assign_g(r,s) \
      (r).c[0]+=(s).c[0]; \
      (r).c[1]+=(s).c[1]; \
      (r).c[2]+=(s).c[2]; \
      (r).c[3]+=(s).c[3]; \
      (r).c[4]+=(s).c[4]; \
      (r).c[5]+=(s).c[5]; \
      (r).c[6]+=(s).c[6]; \
      (r).c[7]+=(s).c[7]; \
      (r).c[8]+=(s).c[8]; \
      (r).c[9]+=(s).c[9]

/* r-=s */
#define _algebra_vector_sub_assign_g(r,s) \
      (r).c[0]-=(s).c[0]; \
      (r).c[1]-=(s).c[1]; \
      (r).c[2]-=(s).c[2]; \
      (r).c[3]-=(s).c[3]; \
      (r).c[4]-=(s).c[4]; \
      (r).c[5]-=(s).c[5]; \
      (r).c[6]-=(s).c[6]; \
      (r).c[7]-=(s).c[7]; \
      (r).c[8]-=(s).c[8]; \
      (r).c[9]-=(s).c[9]

/* r+=k*s (k real) */
#define _algebra_vector_mul_add_assign_g(r,k,s) \
      (r).c[0]+=(k)*(s).c[0]; \
      (r).c[1]+=(k)*(s).c[1]; \
      (r).c[2]+=(k)*(s).c[2]; \
      (r).c[3]+=(k)*(s).c[3]; \
      (r).c[4]+=(k)*(s).c[4]; \
      (r).c[5]+=(k)*(s).c[5]; \
      (r).c[6]+=(k)*(s).c[6]; \
      (r).c[7]+=(k)*(s).c[7]; \
      (r).c[8]+=(k)*(s).c[8]; \
      (r).c[9]+=(k)*(s).c[9]

/* r=k*s (k real) */
#define _algebra_vector_mul_g(r,k,s) \
      (r).c[0]=(k)*(s).c[0]; \
      (r).c[1]=(k)*(s).c[1]; \
      (r).c[2]=(k)*(s).c[2]; \
      (r).c[3]=(k)*(s).c[3]; \
      (r).c[4]=(k)*(s).c[4]; \
      (r).c[5]=(k)*(s).c[5]; \
      (r).c[6]=(k)*(s).c[6]; \
      (r).c[7]=(k)*(s).c[7]; \
      (r).c[8]=(k)*(s).c[8]; \
      (r).c[9]=(k)*(s).c[9]

/* r=0  */
#define _algebra_vector_zero_g(r) \
      (r).c[0]=0.; \
      (r).c[1]=0.; \
      (r).c[2]=0.; \
      (r).c[3]=0.; \
      (r).c[4]=0.; \
      (r).c[5]=0.; \
      (r).c[6]=0.; \
      (r).c[7]=0.; \
      (r).c[8]=0.; \
      (r).c[9]=0.

/* k=|v|^2  */
#define _algebra_vector_sqnorm_g(k,r) \
   (k)=((r).c[0]*(r).c[0])+ \
       ((r).c[1]*(r).c[1])+ \
       ((r).c[2]*(r).c[2])+ \
       ((r).c[3]*(r).c[3])+ \
       ((r).c[4]*(r).c[4])+ \
       ((r).c[5]*(r).c[5])+ \
       ((r).c[6]*(r).c[6])+ \
       ((r).c[7]*(r).c[7])+ \
       ((r).c[8]*(r).c[8])+ \
       ((r).c[9]*(r).c[9])

/* k=Scalar product r*s (r,s algabra vectors)  */
#define _algebra_vector_prod_g(k,r,s) \
   (k)=((r).c[0]*(s).c[0])+ \
       ((r).c[1]*(s).c[1])+ \
       ((r).c[2]*(s).c[2])+ \
       ((r).c[3]*(s).c[3])+ \
       ((r).c[4]*(s).c[4])+ \
       ((r).c[5]*(s).c[5])+ \
       ((r).c[6]*(s).c[6])+ \
       ((r).c[7]*(s).c[7])+ \
       ((r).c[8]*(s).c[8])+ \
       ((r).c[9]*(s).c[9])

/* u = Omega */
#define _symplectic(u) \
   _complex_0((u).c[0]);\
   _complex_0((u).c[1]);\
   _complex_m1((u).c[2]);\
   _complex_0((u).c[3]);\
   _complex_0((u).c[4]);\
   _complex_0((u).c[5]);\
   _complex_0((u).c[6]);\
   _complex_m1((u).c[7])

/* u = u* */
#define _vector_conjugate(u) \
    _complex_star_assign((u).c[0]);\
    _complex_star_assign((u).c[1]);\
    _complex_star_assign((u).c[2]);\
    _complex_star_assign((u).c[3])

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
   _complex_star((u).c[1],(v).c[4]); \
   _complex_minus((u).c[2],(v).c[2]); \
   _complex_minus((u).c[3],(v).c[6]); \
   _complex_star((u).c[4],(v).c[1]); \
   _complex_star((u).c[5],(v).c[5]); \
   _complex_minus((u).c[6],(v).c[3]); \
   _complex_minus((u).c[7],(v).c[7])

/* u=v^T */
#define _suNg_transpose(u,v) \
   (u).c[0]=(v).c[0]; \
   (u).c[1]=(v).c[4]; \
   _complex_star_minus((u).c[2],(v).c[2]); \
   _complex_star_minus((u).c[3],(v).c[6]); \
   (u).c[4]=(v).c[1]; \
   (u).c[5]=(v).c[5]; \
   _complex_star_minus((u).c[6],(v).c[3]); \
   _complex_star_minus((u).c[7],(v).c[7])

/* SP(N) matrix v times SP(N) matrix w */
/* u=v*w */
#define _suNg_times_suNg(u,v,w) \
      _complex_mul((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_assign((u).c[0],(v).c[1],(w).c[4]); \
      _complex_mul_star_massign((u).c[0],(v).c[2],(w).c[2]); \
      _complex_mul_star_massign((u).c[0],(v).c[3],(w).c[6]); \
      _complex_mul((u).c[1],(v).c[0],(w).c[1]);\
      _complex_mul_assign((u).c[1],(v).c[1],(w).c[5]); \
      _complex_mul_star_massign((u).c[1],(v).c[2],(w).c[3]); \
      _complex_mul_star_massign((u).c[1],(v).c[3],(w).c[7]); \
      _complex_mul((u).c[2],(v).c[0],(w).c[2]);\
      _complex_mul_assign((u).c[2],(v).c[1],(w).c[6]); \
      _complex_mul_star_passign((u).c[2],(v).c[2],(w).c[0]); \
      _complex_mul_star_passign((u).c[2],(v).c[3],(w).c[4]); \
      _complex_mul((u).c[3],(v).c[0],(w).c[3]);\
      _complex_mul_assign((u).c[3],(v).c[1],(w).c[7]); \
      _complex_mul_star_passign((u).c[3],(v).c[2],(w).c[1]); \
      _complex_mul_star_passign((u).c[3],(v).c[3],(w).c[5]); \
      _complex_mul((u).c[4],(v).c[4],(w).c[0]);\
      _complex_mul_assign((u).c[4],(v).c[5],(w).c[4]); \
      _complex_mul_star_massign((u).c[4],(v).c[6],(w).c[2]); \
      _complex_mul_star_massign((u).c[4],(v).c[7],(w).c[6]); \
      _complex_mul((u).c[5],(v).c[4],(w).c[1]);\
      _complex_mul_assign((u).c[5],(v).c[5],(w).c[5]); \
      _complex_mul_star_massign((u).c[5],(v).c[6],(w).c[3]); \
      _complex_mul_star_massign((u).c[5],(v).c[7],(w).c[7]); \
      _complex_mul((u).c[6],(v).c[4],(w).c[2]);\
      _complex_mul_assign((u).c[6],(v).c[5],(w).c[6]); \
      _complex_mul_star_passign((u).c[6],(v).c[6],(w).c[0]); \
      _complex_mul_star_passign((u).c[6],(v).c[7],(w).c[4]); \
      _complex_mul((u).c[7],(v).c[4],(w).c[3]);\
      _complex_mul_assign((u).c[7],(v).c[5],(w).c[7]); \
      _complex_mul_star_passign((u).c[7],(v).c[6],(w).c[1]); \
      _complex_mul_star_passign((u).c[7],(v).c[7],(w).c[5]); \
/* u=v*w^+ */
#define _suNg_times_suNg_dagger(u,v,w) \
      _complex_mul_star((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_star_assign((u).c[0],(v).c[1],(w).c[1]); \
      _complex_mul_star_assign((u).c[0],(v).c[2],(w).c[2]); \
      _complex_mul_star_assign((u).c[0],(v).c[3],(w).c[3]); \
      _complex_mul_star((u).c[1],(v).c[0],(w).c[4]);\
      _complex_mul_star_assign((u).c[1],(v).c[1],(w).c[5]); \
      _complex_mul_star_assign((u).c[1],(v).c[2],(w).c[6]); \
      _complex_mul_star_assign((u).c[1],(v).c[3],(w).c[7]); \
      _complex_minus_mul((u).c[2],(v).c[0],(w).c[2]);\
      _complex_mul_massign((u).c[2],(v).c[1],(w).c[3]); \
      _complex_mul_passign((u).c[2],(v).c[2],(w).c[0]); \
      _complex_mul_passign((u).c[2],(v).c[3],(w).c[1]); \
      _complex_minus_mul((u).c[3],(v).c[0],(w).c[6]);\
      _complex_mul_massign((u).c[3],(v).c[1],(w).c[7]); \
      _complex_mul_passign((u).c[3],(v).c[2],(w).c[4]); \
      _complex_mul_passign((u).c[3],(v).c[3],(w).c[5]); \
      _complex_mul_star((u).c[4],(v).c[4],(w).c[0]);\
      _complex_mul_star_assign((u).c[4],(v).c[5],(w).c[1]); \
      _complex_mul_star_assign((u).c[4],(v).c[6],(w).c[2]); \
      _complex_mul_star_assign((u).c[4],(v).c[7],(w).c[3]); \
      _complex_mul_star((u).c[5],(v).c[4],(w).c[4]);\
      _complex_mul_star_assign((u).c[5],(v).c[5],(w).c[5]); \
      _complex_mul_star_assign((u).c[5],(v).c[6],(w).c[6]); \
      _complex_mul_star_assign((u).c[5],(v).c[7],(w).c[7]); \
      _complex_minus_mul((u).c[6],(v).c[4],(w).c[2]);\
      _complex_mul_massign((u).c[6],(v).c[5],(w).c[3]); \
      _complex_mul_passign((u).c[6],(v).c[6],(w).c[0]); \
      _complex_mul_passign((u).c[6],(v).c[7],(w).c[1]); \
      _complex_minus_mul((u).c[7],(v).c[4],(w).c[6]);\
      _complex_mul_massign((u).c[7],(v).c[5],(w).c[7]); \
      _complex_mul_passign((u).c[7],(v).c[6],(w).c[4]); \
      _complex_mul_passign((u).c[7],(v).c[7],(w).c[5]); \
/* u=v^+*w */
#define _suNg_dagger_times_suNg(u,v,w) \
      _complex_mul_star((u).c[0],(w).c[0],(v).c[0]);\
      _complex_mul_star_assign((u).c[0],(w).c[4],(v).c[4]); \
      _complex_mul_star_assign((u).c[0],(v).c[2],(w).c[2]); \
      _complex_mul_star_assign((u).c[0],(v).c[6],(w).c[6]); \
      _complex_mul_star((u).c[1],(w).c[1],(v).c[0]);\
      _complex_mul_star_assign((u).c[1],(w).c[5],(v).c[4]); \
      _complex_mul_star_assign((u).c[1],(v).c[2],(w).c[3]); \
      _complex_mul_star_assign((u).c[1],(v).c[6],(w).c[7]); \
      _complex_mul_star((u).c[2],(w).c[2],(v).c[0]);\
      _complex_mul_star_assign((u).c[2],(w).c[6],(v).c[4]); \
      _complex_mul_star_massign((u).c[2],(v).c[2],(w).c[0]); \
      _complex_mul_star_massign((u).c[2],(v).c[6],(w).c[4]); \
      _complex_mul_star((u).c[3],(w).c[3],(v).c[0]);\
      _complex_mul_star_assign((u).c[3],(w).c[7],(v).c[4]); \
      _complex_mul_star_massign((u).c[3],(v).c[2],(w).c[1]); \
      _complex_mul_star_massign((u).c[3],(v).c[6],(w).c[5]); \
      _complex_mul_star((u).c[4],(w).c[0],(v).c[1]);\
      _complex_mul_star_assign((u).c[4],(w).c[4],(v).c[5]); \
      _complex_mul_star_assign((u).c[4],(v).c[3],(w).c[2]); \
      _complex_mul_star_assign((u).c[4],(v).c[7],(w).c[6]); \
      _complex_mul_star((u).c[5],(w).c[1],(v).c[1]);\
      _complex_mul_star_assign((u).c[5],(w).c[5],(v).c[5]); \
      _complex_mul_star_assign((u).c[5],(v).c[3],(w).c[3]); \
      _complex_mul_star_assign((u).c[5],(v).c[7],(w).c[7]); \
      _complex_mul_star((u).c[6],(w).c[2],(v).c[1]);\
      _complex_mul_star_assign((u).c[6],(w).c[6],(v).c[5]); \
      _complex_mul_star_massign((u).c[6],(v).c[3],(w).c[0]); \
      _complex_mul_star_massign((u).c[6],(v).c[7],(w).c[4]); \
      _complex_mul_star((u).c[7],(w).c[3],(v).c[1]);\
      _complex_mul_star_assign((u).c[7],(w).c[7],(v).c[5]); \
      _complex_mul_star_massign((u).c[7],(v).c[3],(w).c[1]); \
      _complex_mul_star_massign((u).c[7],(v).c[7],(w).c[5])

/* u=0 */
#define _suNg_zero(u) \
    _complex_0((u).c[0]);\
    _complex_0((u).c[1]);\
    _complex_0((u).c[2]);\
    _complex_0((u).c[3]);\
    _complex_0((u).c[4]);\
    _complex_0((u).c[5]);\
    _complex_0((u).c[6]);\
    _complex_0((u).c[7])

/* u=1 */
#define _suNg_unit(u) \
   _complex_1((u).c[0]);\
   _complex_0((u).c[1]);\
   _complex_0((u).c[2]);\
   _complex_0((u).c[3]);\
   _complex_0((u).c[4]);\
   _complex_1((u).c[5]);\
   _complex_0((u).c[6]);\
   _complex_0((u).c[7]);\

/* Expand a symplectic matrix into a full matrix */
#define _suNg_expand( v, u ) \
   (v).c[0].re = (u).c[0].re;\
   (v).c[0].im = (u).c[0].im;\
   (v).c[10].re = (u).c[0].re;\
   (v).c[10].im =-(u).c[0].im;\
   (v).c[2].re = (u).c[2].re;\
   (v).c[2].im = (u).c[2].im;\
   (v).c[8].re =-(u).c[2].re;\
   (v).c[8].im = (u).c[2].im;\
   (v).c[1].re = (u).c[1].re;\
   (v).c[1].im = (u).c[1].im;\
   (v).c[11].re = (u).c[1].re;\
   (v).c[11].im =-(u).c[1].im;\
   (v).c[3].re = (u).c[3].re;\
   (v).c[3].im = (u).c[3].im;\
   (v).c[9].re =-(u).c[3].re;\
   (v).c[9].im = (u).c[3].im;\
   (v).c[4].re = (u).c[4].re;\
   (v).c[4].im = (u).c[4].im;\
   (v).c[14].re = (u).c[4].re;\
   (v).c[14].im =-(u).c[4].im;\
   (v).c[6].re = (u).c[6].re;\
   (v).c[6].im = (u).c[6].im;\
   (v).c[12].re =-(u).c[6].re;\
   (v).c[12].im = (u).c[6].im;\
   (v).c[5].re = (u).c[5].re;\
   (v).c[5].im = (u).c[5].im;\
   (v).c[15].re = (u).c[5].re;\
   (v).c[15].im =-(u).c[5].im;\
   (v).c[7].re = (u).c[7].re;\
   (v).c[7].im = (u).c[7].im;\
   (v).c[13].re =-(u).c[7].re;\
   (v).c[13].im = (u).c[7].im;\

/* u=-v */
#define _suNg_minus(u,v) \
   _complex_minus((u).c[0],(v).c[0]);\
   _complex_minus((u).c[1],(v).c[1]);\
   _complex_minus((u).c[2],(v).c[2]);\
   _complex_minus((u).c[3],(v).c[3]);\
   _complex_minus((u).c[4],(v).c[4]);\
   _complex_minus((u).c[5],(v).c[5]);\
   _complex_minus((u).c[6],(v).c[6]);\
   _complex_minus((u).c[7],(v).c[7])

/* u=r*v (u,v mat, r real) */
#define _suNg_mul(u,r,v) \
   _complex_mulr((u).c[0],(r),(v).c[0]);\
   _complex_mulr((u).c[1],(r),(v).c[1]);\
   _complex_mulr((u).c[2],(r),(v).c[2]);\
   _complex_mulr((u).c[3],(r),(v).c[3]);\
   _complex_mulr((u).c[4],(r),(v).c[4]);\
   _complex_mulr((u).c[5],(r),(v).c[5]);\
   _complex_mulr((u).c[6],(r),(v).c[6]);\
   _complex_mulr((u).c[7],(r),(v).c[7])

/* u=r*v (u,v mat, r complex) */
#define _suNg_mulc(u,r,v) \
   _complex_mul((u).c[0],(r),(v).c[0]);\
   _complex_mul((u).c[1],(r),(v).c[1]);\
   _complex_mul((u).c[2],(r),(v).c[2]);\
   _complex_mul((u).c[3],(r),(v).c[3]);\
   _complex_mul((u).c[4],(r),(v).c[4]);\
   _complex_mul((u).c[5],(r),(v).c[5]);\
   _complex_mul((u).c[6],(r),(v).c[6]);\
   _complex_mul((u).c[7],(r),(v).c[7])

/* u+=v */
#define _suNg_add_assign(u,v) \
   _complex_add_assign((u).c[0],(v).c[0]);\
   _complex_add_assign((u).c[1],(v).c[1]);\
   _complex_add_assign((u).c[2],(v).c[2]);\
   _complex_add_assign((u).c[3],(v).c[3]);\
   _complex_add_assign((u).c[4],(v).c[4]);\
   _complex_add_assign((u).c[5],(v).c[5]);\
   _complex_add_assign((u).c[6],(v).c[6]);\
   _complex_add_assign((u).c[7],(v).c[7])

/* u-=v */
#define _suNg_sub_assign(u,v) \
   _complex_sub_assign((u).c[0],(v).c[0]);\
   _complex_sub_assign((u).c[1],(v).c[1]);\
   _complex_sub_assign((u).c[2],(v).c[2]);\
   _complex_sub_assign((u).c[3],(v).c[3]);\
   _complex_sub_assign((u).c[4],(v).c[4]);\
   _complex_sub_assign((u).c[5],(v).c[5]);\
   _complex_sub_assign((u).c[6],(v).c[6]);\
   _complex_sub_assign((u).c[7],(v).c[7])

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
   (k)+=_complex_prod_re((u).c[7],(u).c[7]);\
      (k)*=2;\
/* k=| 1 - u |2 ) */
#define _suNg_sqnorm_m1(k,u) \
   (k)=\
    +_complex_prod_m1_re((u).c[0],(u).c[0])\
    +_complex_prod_re((u).c[1],(u).c[1])\
    +_complex_prod_re((u).c[2],(u).c[2])\
    +_complex_prod_re((u).c[3],(u).c[3])\
    +_complex_prod_re((u).c[4],(u).c[4])\
    +_complex_prod_m1_re((u).c[5],(u).c[5])\
    +_complex_prod_re((u).c[6],(u).c[6])\
    +_complex_prod_re((u).c[7],(u).c[7]);\
      (k)*=2;

/* k=Re Tr (u) */
#define _suNg_trace_re(k,u) \
   (k)=_complex_re((u).c[0])+ \
       _complex_re((u).c[5]);\
      (k)*=2;

/* k=Im Tr (u) */
#define _suNg_trace_im(k,u) \
   (k)=0;

/* This fuction computes the hmc force matrix */
/* this fuction accumulates on the original matrix u */
#define _suNg_FMAT(u,s) \
   _complex_mul_star_assign((u).c[0],(s).c[0].c[0],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[0],(s).c[1].c[0],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[1],(s).c[0].c[0],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[1],(s).c[1].c[0],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[2],(s).c[0].c[0],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[2],(s).c[1].c[0],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[3],(s).c[0].c[0],(s).c[2].c[3]); \
   _complex_mul_star_assign((u).c[3],(s).c[1].c[0],(s).c[3].c[3]);\
   _complex_mul_star_assign((u).c[4],(s).c[0].c[1],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[4],(s).c[1].c[1],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[5],(s).c[0].c[1],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[5],(s).c[1].c[1],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[6],(s).c[0].c[1],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[6],(s).c[1].c[1],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[7],(s).c[0].c[1],(s).c[2].c[3]); \
   _complex_mul_star_assign((u).c[7],(s).c[1].c[1],(s).c[3].c[3]);\
   _complex_mul_star_assign((u).c[8],(s).c[0].c[2],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[8],(s).c[1].c[2],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[9],(s).c[0].c[2],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[9],(s).c[1].c[2],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[10],(s).c[0].c[2],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[10],(s).c[1].c[2],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[11],(s).c[0].c[2],(s).c[2].c[3]); \
   _complex_mul_star_assign((u).c[11],(s).c[1].c[2],(s).c[3].c[3]);\
   _complex_mul_star_assign((u).c[12],(s).c[0].c[3],(s).c[2].c[0]); \
   _complex_mul_star_assign((u).c[12],(s).c[1].c[3],(s).c[3].c[0]);\
   _complex_mul_star_assign((u).c[13],(s).c[0].c[3],(s).c[2].c[1]); \
   _complex_mul_star_assign((u).c[13],(s).c[1].c[3],(s).c[3].c[1]);\
   _complex_mul_star_assign((u).c[14],(s).c[0].c[3],(s).c[2].c[2]); \
   _complex_mul_star_assign((u).c[14],(s).c[1].c[3],(s).c[3].c[2]);\
   _complex_mul_star_assign((u).c[15],(s).c[0].c[3],(s).c[2].c[3]); \
   _complex_mul_star_assign((u).c[15],(s).c[1].c[3],(s).c[3].c[3])

/* u=v*w */
#define _suNgfull_times_suNgfull(u,v,w) \
      _complex_mul((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_assign((u).c[0],(v).c[1],(w).c[4]); \
      _complex_mul_assign((u).c[0],(v).c[2],(w).c[8]); \
      _complex_mul_assign((u).c[0],(v).c[3],(w).c[12]); \
      _complex_mul((u).c[1],(v).c[0],(w).c[1]);\
      _complex_mul_assign((u).c[1],(v).c[1],(w).c[5]); \
      _complex_mul_assign((u).c[1],(v).c[2],(w).c[9]); \
      _complex_mul_assign((u).c[1],(v).c[3],(w).c[13]); \
      _complex_mul((u).c[2],(v).c[0],(w).c[2]);\
      _complex_mul_assign((u).c[2],(v).c[1],(w).c[6]); \
      _complex_mul_assign((u).c[2],(v).c[2],(w).c[10]); \
      _complex_mul_assign((u).c[2],(v).c[3],(w).c[14]); \
      _complex_mul((u).c[3],(v).c[0],(w).c[3]);\
      _complex_mul_assign((u).c[3],(v).c[1],(w).c[7]); \
      _complex_mul_assign((u).c[3],(v).c[2],(w).c[11]); \
      _complex_mul_assign((u).c[3],(v).c[3],(w).c[15]); \
      _complex_mul((u).c[4],(v).c[4],(w).c[0]);\
      _complex_mul_assign((u).c[4],(v).c[5],(w).c[4]); \
      _complex_mul_assign((u).c[4],(v).c[6],(w).c[8]); \
      _complex_mul_assign((u).c[4],(v).c[7],(w).c[12]); \
      _complex_mul((u).c[5],(v).c[4],(w).c[1]);\
      _complex_mul_assign((u).c[5],(v).c[5],(w).c[5]); \
      _complex_mul_assign((u).c[5],(v).c[6],(w).c[9]); \
      _complex_mul_assign((u).c[5],(v).c[7],(w).c[13]); \
      _complex_mul((u).c[6],(v).c[4],(w).c[2]);\
      _complex_mul_assign((u).c[6],(v).c[5],(w).c[6]); \
      _complex_mul_assign((u).c[6],(v).c[6],(w).c[10]); \
      _complex_mul_assign((u).c[6],(v).c[7],(w).c[14]); \
      _complex_mul((u).c[7],(v).c[4],(w).c[3]);\
      _complex_mul_assign((u).c[7],(v).c[5],(w).c[7]); \
      _complex_mul_assign((u).c[7],(v).c[6],(w).c[11]); \
      _complex_mul_assign((u).c[7],(v).c[7],(w).c[15]); \
      _complex_mul((u).c[8],(v).c[8],(w).c[0]);\
      _complex_mul_assign((u).c[8],(v).c[9],(w).c[4]); \
      _complex_mul_assign((u).c[8],(v).c[10],(w).c[8]); \
      _complex_mul_assign((u).c[8],(v).c[11],(w).c[12]); \
      _complex_mul((u).c[9],(v).c[8],(w).c[1]);\
      _complex_mul_assign((u).c[9],(v).c[9],(w).c[5]); \
      _complex_mul_assign((u).c[9],(v).c[10],(w).c[9]); \
      _complex_mul_assign((u).c[9],(v).c[11],(w).c[13]); \
      _complex_mul((u).c[10],(v).c[8],(w).c[2]);\
      _complex_mul_assign((u).c[10],(v).c[9],(w).c[6]); \
      _complex_mul_assign((u).c[10],(v).c[10],(w).c[10]); \
      _complex_mul_assign((u).c[10],(v).c[11],(w).c[14]); \
      _complex_mul((u).c[11],(v).c[8],(w).c[3]);\
      _complex_mul_assign((u).c[11],(v).c[9],(w).c[7]); \
      _complex_mul_assign((u).c[11],(v).c[10],(w).c[11]); \
      _complex_mul_assign((u).c[11],(v).c[11],(w).c[15]); \
      _complex_mul((u).c[12],(v).c[12],(w).c[0]);\
      _complex_mul_assign((u).c[12],(v).c[13],(w).c[4]); \
      _complex_mul_assign((u).c[12],(v).c[14],(w).c[8]); \
      _complex_mul_assign((u).c[12],(v).c[15],(w).c[12]); \
      _complex_mul((u).c[13],(v).c[12],(w).c[1]);\
      _complex_mul_assign((u).c[13],(v).c[13],(w).c[5]); \
      _complex_mul_assign((u).c[13],(v).c[14],(w).c[9]); \
      _complex_mul_assign((u).c[13],(v).c[15],(w).c[13]); \
      _complex_mul((u).c[14],(v).c[12],(w).c[2]);\
      _complex_mul_assign((u).c[14],(v).c[13],(w).c[6]); \
      _complex_mul_assign((u).c[14],(v).c[14],(w).c[10]); \
      _complex_mul_assign((u).c[14],(v).c[15],(w).c[14]); \
      _complex_mul((u).c[15],(v).c[12],(w).c[3]);\
      _complex_mul_assign((u).c[15],(v).c[13],(w).c[7]); \
      _complex_mul_assign((u).c[15],(v).c[14],(w).c[11]); \
      _complex_mul_assign((u).c[15],(v).c[15],(w).c[15])

/* u=v*w^+ */
#define _suNgfull_times_suNgfull_dagger(u,v,w) \
      _complex_mul_star((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_star_assign((u).c[0],(v).c[1],(w).c[1]); \
      _complex_mul_star_assign((u).c[0],(v).c[2],(w).c[2]); \
      _complex_mul_star_assign((u).c[0],(v).c[3],(w).c[3]); \
      _complex_mul_star((u).c[1],(v).c[0],(w).c[4]);\
      _complex_mul_star_assign((u).c[1],(v).c[1],(w).c[5]); \
      _complex_mul_star_assign((u).c[1],(v).c[2],(w).c[6]); \
      _complex_mul_star_assign((u).c[1],(v).c[3],(w).c[7]); \
      _complex_mul_star((u).c[2],(v).c[0],(w).c[8]);\
      _complex_mul_star_assign((u).c[2],(v).c[1],(w).c[9]); \
      _complex_mul_star_assign((u).c[2],(v).c[2],(w).c[10]); \
      _complex_mul_star_assign((u).c[2],(v).c[3],(w).c[11]); \
      _complex_mul_star((u).c[3],(v).c[0],(w).c[12]);\
      _complex_mul_star_assign((u).c[3],(v).c[1],(w).c[13]); \
      _complex_mul_star_assign((u).c[3],(v).c[2],(w).c[14]); \
      _complex_mul_star_assign((u).c[3],(v).c[3],(w).c[15]); \
      _complex_mul_star((u).c[4],(v).c[4],(w).c[0]);\
      _complex_mul_star_assign((u).c[4],(v).c[5],(w).c[1]); \
      _complex_mul_star_assign((u).c[4],(v).c[6],(w).c[2]); \
      _complex_mul_star_assign((u).c[4],(v).c[7],(w).c[3]); \
      _complex_mul_star((u).c[5],(v).c[4],(w).c[4]);\
      _complex_mul_star_assign((u).c[5],(v).c[5],(w).c[5]); \
      _complex_mul_star_assign((u).c[5],(v).c[6],(w).c[6]); \
      _complex_mul_star_assign((u).c[5],(v).c[7],(w).c[7]); \
      _complex_mul_star((u).c[6],(v).c[4],(w).c[8]);\
      _complex_mul_star_assign((u).c[6],(v).c[5],(w).c[9]); \
      _complex_mul_star_assign((u).c[6],(v).c[6],(w).c[10]); \
      _complex_mul_star_assign((u).c[6],(v).c[7],(w).c[11]); \
      _complex_mul_star((u).c[7],(v).c[4],(w).c[12]);\
      _complex_mul_star_assign((u).c[7],(v).c[5],(w).c[13]); \
      _complex_mul_star_assign((u).c[7],(v).c[6],(w).c[14]); \
      _complex_mul_star_assign((u).c[7],(v).c[7],(w).c[15]); \
      _complex_mul_star((u).c[8],(v).c[8],(w).c[0]);\
      _complex_mul_star_assign((u).c[8],(v).c[9],(w).c[1]); \
      _complex_mul_star_assign((u).c[8],(v).c[10],(w).c[2]); \
      _complex_mul_star_assign((u).c[8],(v).c[11],(w).c[3]); \
      _complex_mul_star((u).c[9],(v).c[8],(w).c[4]);\
      _complex_mul_star_assign((u).c[9],(v).c[9],(w).c[5]); \
      _complex_mul_star_assign((u).c[9],(v).c[10],(w).c[6]); \
      _complex_mul_star_assign((u).c[9],(v).c[11],(w).c[7]); \
      _complex_mul_star((u).c[10],(v).c[8],(w).c[8]);\
      _complex_mul_star_assign((u).c[10],(v).c[9],(w).c[9]); \
      _complex_mul_star_assign((u).c[10],(v).c[10],(w).c[10]); \
      _complex_mul_star_assign((u).c[10],(v).c[11],(w).c[11]); \
      _complex_mul_star((u).c[11],(v).c[8],(w).c[12]);\
      _complex_mul_star_assign((u).c[11],(v).c[9],(w).c[13]); \
      _complex_mul_star_assign((u).c[11],(v).c[10],(w).c[14]); \
      _complex_mul_star_assign((u).c[11],(v).c[11],(w).c[15]); \
      _complex_mul_star((u).c[12],(v).c[12],(w).c[0]);\
      _complex_mul_star_assign((u).c[12],(v).c[13],(w).c[1]); \
      _complex_mul_star_assign((u).c[12],(v).c[14],(w).c[2]); \
      _complex_mul_star_assign((u).c[12],(v).c[15],(w).c[3]); \
      _complex_mul_star((u).c[13],(v).c[12],(w).c[4]);\
      _complex_mul_star_assign((u).c[13],(v).c[13],(w).c[5]); \
      _complex_mul_star_assign((u).c[13],(v).c[14],(w).c[6]); \
      _complex_mul_star_assign((u).c[13],(v).c[15],(w).c[7]); \
      _complex_mul_star((u).c[14],(v).c[12],(w).c[8]);\
      _complex_mul_star_assign((u).c[14],(v).c[13],(w).c[9]); \
      _complex_mul_star_assign((u).c[14],(v).c[14],(w).c[10]); \
      _complex_mul_star_assign((u).c[14],(v).c[15],(w).c[11]); \
      _complex_mul_star((u).c[15],(v).c[12],(w).c[12]);\
      _complex_mul_star_assign((u).c[15],(v).c[13],(w).c[13]); \
      _complex_mul_star_assign((u).c[15],(v).c[14],(w).c[14]); \
      _complex_mul_star_assign((u).c[15],(v).c[15],(w).c[15])

/* u=v^+*w */
#define _suNgfull_dagger_times_suNgfull(u,v,w) \
      _complex_mul_star((u).c[0],(w).c[0],(v).c[0]);\
      _complex_mul_star_assign((u).c[0],(w).c[4],(v).c[4]); \
      _complex_mul_star_assign((u).c[0],(w).c[8],(v).c[8]); \
      _complex_mul_star_assign((u).c[0],(w).c[12],(v).c[12]); \
      _complex_mul_star((u).c[1],(w).c[1],(v).c[0]);\
      _complex_mul_star_assign((u).c[1],(w).c[5],(v).c[4]); \
      _complex_mul_star_assign((u).c[1],(w).c[9],(v).c[8]); \
      _complex_mul_star_assign((u).c[1],(w).c[13],(v).c[12]); \
      _complex_mul_star((u).c[2],(w).c[2],(v).c[0]);\
      _complex_mul_star_assign((u).c[2],(w).c[6],(v).c[4]); \
      _complex_mul_star_assign((u).c[2],(w).c[10],(v).c[8]); \
      _complex_mul_star_assign((u).c[2],(w).c[14],(v).c[12]); \
      _complex_mul_star((u).c[3],(w).c[3],(v).c[0]);\
      _complex_mul_star_assign((u).c[3],(w).c[7],(v).c[4]); \
      _complex_mul_star_assign((u).c[3],(w).c[11],(v).c[8]); \
      _complex_mul_star_assign((u).c[3],(w).c[15],(v).c[12]); \
      _complex_mul_star((u).c[4],(w).c[0],(v).c[1]);\
      _complex_mul_star_assign((u).c[4],(w).c[4],(v).c[5]); \
      _complex_mul_star_assign((u).c[4],(w).c[8],(v).c[9]); \
      _complex_mul_star_assign((u).c[4],(w).c[12],(v).c[13]); \
      _complex_mul_star((u).c[5],(w).c[1],(v).c[1]);\
      _complex_mul_star_assign((u).c[5],(w).c[5],(v).c[5]); \
      _complex_mul_star_assign((u).c[5],(w).c[9],(v).c[9]); \
      _complex_mul_star_assign((u).c[5],(w).c[13],(v).c[13]); \
      _complex_mul_star((u).c[6],(w).c[2],(v).c[1]);\
      _complex_mul_star_assign((u).c[6],(w).c[6],(v).c[5]); \
      _complex_mul_star_assign((u).c[6],(w).c[10],(v).c[9]); \
      _complex_mul_star_assign((u).c[6],(w).c[14],(v).c[13]); \
      _complex_mul_star((u).c[7],(w).c[3],(v).c[1]);\
      _complex_mul_star_assign((u).c[7],(w).c[7],(v).c[5]); \
      _complex_mul_star_assign((u).c[7],(w).c[11],(v).c[9]); \
      _complex_mul_star_assign((u).c[7],(w).c[15],(v).c[13]); \
      _complex_mul_star((u).c[8],(w).c[0],(v).c[2]);\
      _complex_mul_star_assign((u).c[8],(w).c[4],(v).c[6]); \
      _complex_mul_star_assign((u).c[8],(w).c[8],(v).c[10]); \
      _complex_mul_star_assign((u).c[8],(w).c[12],(v).c[14]); \
      _complex_mul_star((u).c[9],(w).c[1],(v).c[2]);\
      _complex_mul_star_assign((u).c[9],(w).c[5],(v).c[6]); \
      _complex_mul_star_assign((u).c[9],(w).c[9],(v).c[10]); \
      _complex_mul_star_assign((u).c[9],(w).c[13],(v).c[14]); \
      _complex_mul_star((u).c[10],(w).c[2],(v).c[2]);\
      _complex_mul_star_assign((u).c[10],(w).c[6],(v).c[6]); \
      _complex_mul_star_assign((u).c[10],(w).c[10],(v).c[10]); \
      _complex_mul_star_assign((u).c[10],(w).c[14],(v).c[14]); \
      _complex_mul_star((u).c[11],(w).c[3],(v).c[2]);\
      _complex_mul_star_assign((u).c[11],(w).c[7],(v).c[6]); \
      _complex_mul_star_assign((u).c[11],(w).c[11],(v).c[10]); \
      _complex_mul_star_assign((u).c[11],(w).c[15],(v).c[14]); \
      _complex_mul_star((u).c[12],(w).c[0],(v).c[3]);\
      _complex_mul_star_assign((u).c[12],(w).c[4],(v).c[7]); \
      _complex_mul_star_assign((u).c[12],(w).c[8],(v).c[11]); \
      _complex_mul_star_assign((u).c[12],(w).c[12],(v).c[15]); \
      _complex_mul_star((u).c[13],(w).c[1],(v).c[3]);\
      _complex_mul_star_assign((u).c[13],(w).c[5],(v).c[7]); \
      _complex_mul_star_assign((u).c[13],(w).c[9],(v).c[11]); \
      _complex_mul_star_assign((u).c[13],(w).c[13],(v).c[15]); \
      _complex_mul_star((u).c[14],(w).c[2],(v).c[3]);\
      _complex_mul_star_assign((u).c[14],(w).c[6],(v).c[7]); \
      _complex_mul_star_assign((u).c[14],(w).c[10],(v).c[11]); \
      _complex_mul_star_assign((u).c[14],(w).c[14],(v).c[15]); \
      _complex_mul_star((u).c[15],(w).c[3],(v).c[3]);\
      _complex_mul_star_assign((u).c[15],(w).c[7],(v).c[7]); \
      _complex_mul_star_assign((u).c[15],(w).c[11],(v).c[11]); \
      _complex_mul_star_assign((u).c[15],(w).c[15],(v).c[15])

/* u+=v */
#define _suNgfull_add_assign(u,v) \
   _complex_add_assign((u).c[0],(v).c[0]);\
   _complex_add_assign((u).c[1],(v).c[1]);\
   _complex_add_assign((u).c[2],(v).c[2]);\
   _complex_add_assign((u).c[3],(v).c[3]);\
   _complex_add_assign((u).c[4],(v).c[4]);\
   _complex_add_assign((u).c[5],(v).c[5]);\
   _complex_add_assign((u).c[6],(v).c[6]);\
   _complex_add_assign((u).c[7],(v).c[7]);\
   _complex_add_assign((u).c[8],(v).c[8]);\
   _complex_add_assign((u).c[9],(v).c[9]);\
   _complex_add_assign((u).c[10],(v).c[10]);\
   _complex_add_assign((u).c[11],(v).c[11]);\
   _complex_add_assign((u).c[12],(v).c[12]);\
   _complex_add_assign((u).c[13],(v).c[13]);\
   _complex_add_assign((u).c[14],(v).c[14]);\
   _complex_add_assign((u).c[15],(v).c[15])

/* u-=v */
#define _suNgfull_sub_assign(u,v) \
   _complex_sub_assign((u).c[0],(v).c[0]);\
   _complex_sub_assign((u).c[1],(v).c[1]);\
   _complex_sub_assign((u).c[2],(v).c[2]);\
   _complex_sub_assign((u).c[3],(v).c[3]);\
   _complex_sub_assign((u).c[4],(v).c[4]);\
   _complex_sub_assign((u).c[5],(v).c[5]);\
   _complex_sub_assign((u).c[6],(v).c[6]);\
   _complex_sub_assign((u).c[7],(v).c[7]);\
   _complex_sub_assign((u).c[8],(v).c[8]);\
   _complex_sub_assign((u).c[9],(v).c[9]);\
   _complex_sub_assign((u).c[10],(v).c[10]);\
   _complex_sub_assign((u).c[11],(v).c[11]);\
   _complex_sub_assign((u).c[12],(v).c[12]);\
   _complex_sub_assign((u).c[13],(v).c[13]);\
   _complex_sub_assign((u).c[14],(v).c[14]);\
   _complex_sub_assign((u).c[15],(v).c[15])

/* u=0 */
#define _suNgfull_zero(u) \
    _complex_0((u).c[0]);\
    _complex_0((u).c[1]);\
    _complex_0((u).c[2]);\
    _complex_0((u).c[3]);\
    _complex_0((u).c[4]);\
    _complex_0((u).c[5]);\
    _complex_0((u).c[6]);\
    _complex_0((u).c[7]);\
    _complex_0((u).c[8]);\
    _complex_0((u).c[9]);\
    _complex_0((u).c[10]);\
    _complex_0((u).c[11]);\
    _complex_0((u).c[12]);\
    _complex_0((u).c[13]);\
    _complex_0((u).c[14]);\
    _complex_0((u).c[15])

/* u=1 */
#define _suNgfull_unit(u) \
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

/* u=r*v (u,v mat, r real) */
#define _suNgfull_mul(u,r,v) \
   _complex_mulr((u).c[0],(r),(v).c[0]);\
   _complex_mulr((u).c[1],(r),(v).c[1]);\
   _complex_mulr((u).c[2],(r),(v).c[2]);\
   _complex_mulr((u).c[3],(r),(v).c[3]);\
   _complex_mulr((u).c[4],(r),(v).c[4]);\
   _complex_mulr((u).c[5],(r),(v).c[5]);\
   _complex_mulr((u).c[6],(r),(v).c[6]);\
   _complex_mulr((u).c[7],(r),(v).c[7]);\
   _complex_mulr((u).c[8],(r),(v).c[8]);\
   _complex_mulr((u).c[9],(r),(v).c[9]);\
   _complex_mulr((u).c[10],(r),(v).c[10]);\
   _complex_mulr((u).c[11],(r),(v).c[11]);\
   _complex_mulr((u).c[12],(r),(v).c[12]);\
   _complex_mulr((u).c[13],(r),(v).c[13]);\
   _complex_mulr((u).c[14],(r),(v).c[14]);\
   _complex_mulr((u).c[15],(r),(v).c[15])

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
    _complex_0((u).c[8]);\
    _complex_0((u).c[9]);\
    _complex_0((u).c[10]);\
    _complex_0((u).c[11]);\
    _complex_0((u).c[12]);\
    _complex_0((u).c[13]);\
    _complex_0((u).c[14]);\
    _complex_0((u).c[15])

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
      int __iz=(iy)+((x)*4)*(stride); \
      (v).c[0]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex_flt*)(in))[__iz]; \
   } while (0) 

#define _suNg_read_spinor_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*4)*(stride); \
      (v).c[0]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex*)(in))[__iz]; \
   } while (0) 

/* Write spinor field component to GPU memory */
/* (input) v = suNg_vector ; (output) out = suNg_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNg_write_spinor_flt_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*4)*(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[3]; \
   } while (0) 

#define _suNg_write_spinor_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*4)*(stride); \
      ((complex*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[3]; \
   } while (0) 

/* Read an suN algebra vector from GPU memory */
/* (output) v = suN_algebra_vector ; (input) in = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_av_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      (v).c[14]=((float*)(in))[__iz]; \
   } while (0) 

#define _suNg_av_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      (v).c[14]=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN algebra vector to GPU memory */
/* (input) v = suN_algebra_vector ; (output) out = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_av_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      ((float*)(out))[__iz]=(v).c[14]; \
   } while (0) 

#define _suNg_av_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      ((double*)(out))[__iz]=(v).c[14]; \
   } while (0) 

/* Mul_add_assign on a suN algebra vector on GPU  */
/* (in/out) v = suN_algebra_vector* ; (input) in = suN_algebra_vector */
/* (input) iy = site ; (input) x = 0..3 direction; (input) r = real */
#define _algebra_vector_mul_add_assign_gpu_g_flt(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      ((float*)(v))[__iz]+=(in).c[14]*(r); \
   } while (0) 

#define _algebra_vector_mul_add_assign_gpu_g(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*15)*(stride); \
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
      ((double*)(v))[__iz]+=(in).c[14]*(r); \
   } while (0) 

/* Read an suN matrix from GPU memory */
/* (output) v = suN ; (input) in = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*32)*(stride); \
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
      (v).c[15].im=((float*)(in))[__iz]; \
   } while (0) 

#define _suNg_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*32)*(stride); \
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
      (v).c[15].im=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN matrix to GPU memory */
/* (input) v = suN ; (output) out = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNg_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*32)*(stride); \
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
      ((float*)(out))[__iz]=(v).c[15].im; \
   } while (0) 

#define _suNg_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*32)*(stride); \
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
      ((double*)(out))[__iz]=(v).c[15].im; \
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
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_1((r).c[_i]); ++_i;\
         _complex_1((r).c[_i]); ++_i;\
         _complex_1((r).c[_i]); ++_i;\
         _complex_1((r).c[_i]); ++_i;\
      }\
      _complex_1((r).c[_i]); ++_i;\
      _complex_1((r).c[_i]); ++_i;\
   } while(0) 

/* r=0 */
#define _vector_zero_f(r) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_0((r).c[_i]); ++_i;\
         _complex_0((r).c[_i]); ++_i;\
         _complex_0((r).c[_i]); ++_i;\
         _complex_0((r).c[_i]); ++_i;\
      }\
      _complex_0((r).c[_i]); ++_i;\
      _complex_0((r).c[_i]); ++_i;\
   } while(0) 

/* r=-s */
#define _vector_minus_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
      }\
      _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
      _complex_minus((r).c[_i],(s).c[_i]); ++_i; \
   } while(0) 

/* r= i*s */
#define _vector_i_plus_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
      }\
      _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
      _complex_i_plus((r).c[_i],(s).c[_i]); ++_i; \
   } while(0) 

/* r=-i*s */
#define _vector_i_minus_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
         _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
      }\
      _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
      _complex_i_minus((r).c[_i],(s).c[_i]); ++_i; \
   } while(0) 

/* r=k*s (k real) */
#define _vector_mul_f(r,k,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
      }\
      _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
      _complex_mulr((r).c[_i],(k),(s).c[_i]); ++_i;\
   } while(0) 

/* r=z*s (z complex) */
#define _vector_mulc_f(r,z,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
      }\
      _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
      _complex_mul((r).c[_i],(z),(s).c[_i]); ++_i;\
   } while(0) 

/* r=(z^+)*s (z complex) */
#define _vector_mulc_star_f(r,z,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
         _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
         _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
         _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
      }\
      _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
      _complex_mul_star((r).c[_i],(s).c[_i],(z)); ++_i;\
   } while(0) 

/* r=s1+s2 */
#define _vector_add_f(r,s1,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      }\
      _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      _complex_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
   } while(0) 

/* r=s1-s2 */
#define _vector_sub_f(r,s1,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      }\
      _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      _complex_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
   } while(0) 

/* r=s1+i*s2 */
#define _vector_i_add_f(r,s1,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      }\
      _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      _complex_i_add((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
   } while(0) 

/* r=s1-i*s2 */
#define _vector_i_sub_f(r,s1,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
         _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      }\
      _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
      _complex_i_sub((r).c[_i],(s1).c[_i],(s2).c[_i]); ++_i;\
   } while(0) 

/* r+=s */
#define _vector_add_assign_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
      }\
      _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
      _complex_add_assign((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* r-=s */
#define _vector_sub_assign_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
      }\
      _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
      _complex_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* r+=i*s */
#define _vector_i_add_assign_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
      }\
      _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
      _complex_i_add_assign((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* r-=i*s */
#define _vector_i_sub_assign_f(r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
         _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
      }\
      _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
      _complex_i_sub_assign((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k=Re(r^*s) */
#define _vector_prod_re_f(k,r,s) \
   do { \
      int _i;\
      (k)=_complex_prod_re((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k=Im(r*s) */
#define _vector_prod_im_f(k,r,s) \
   do { \
      int _i;\
      (k)=_complex_prod_im((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* r+=z*s (z complex) */
#define _vector_mulc_add_assign_f(r,z,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
      }\
      _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
      _complex_mul_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
   } while(0) 

/* r+=k*s (k real) */
#define _vector_mul_add_assign_f(r,k,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
         _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
      }\
      _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
      _complex_mulr_assign((r).c[_i],(k),(s).c[_i]); ++_i;\
   } while(0) 

/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_f(r,k1,s1,k2,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
      }\
      _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
      _complex_rlc((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
   } while(0) 

/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */
#define _vector_lc_add_assign_f(r,k1,s1,k2,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
         _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
      }\
      _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
      _complex_rlc_assign((r).c[_i],(k1),(s1).c[_i],(k2),(s2).c[_i]); ++_i;\
   } while(0) 

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_f(r,z1,s1,z2,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
      }\
      _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
      _complex_clc((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
   } while(0) 

/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */
#define _vector_clc_add_assign_f(r,z1,s1,z2,s2) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
         _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
      }\
      _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
      _complex_clc_assign((r).c[_i],(z1),(s1).c[_i],(z2),(s2).c[_i]); ++_i;\
   } while(0) 

/* z+=r^*s (c complex) */
#define _vector_prod_assign_f(z,r,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
         _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
         _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
         _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
      }\
      _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
      _complex_prod_assign((z),(r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k+=Re(r^*s) */
#define _vector_prod_add_assign_re_f(k,r,s) \
   do { \
      int _i;\
      (k)+=_complex_prod_re((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)+=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k+=Im(r*s) */
#define _vector_prod_add_assign_im_f(k,r,s) \
   do { \
      int _i;\
      (k)+=_complex_prod_im((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)+=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k-=Re(r^*s) */
#define _vector_prod_sub_assign_re_f(k,r,s) \
   do { \
      int _i;\
      (k)-=_complex_prod_re((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)-=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)-=_complex_prod_re((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* k-=Im(r*s) */
#define _vector_prod_sub_assign_im_f(k,r,s) \
   do { \
      int _i;\
      (k)-=_complex_prod_im((r).c[0],(s).c[0]);\
      for (_i=1; _i<8; ){\
         (k)-=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
         (k)-=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
      }\
      (k)-=_complex_prod_im((r).c[_i],(s).c[_i]); ++_i;\
   } while(0) 

/* r-=z*s (z complex) */
#define _vector_project_f(r,z,s) \
   do { \
      int _i;for (_i=0; _i<8; ){\
         _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
         _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
      }\
      _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
      _complex_mul_sub_assign((r).c[_i],(z),(s).c[_i]); ++_i;\
   } while(0) 

/* SU(N) matrix u times SU(N) vector s */
/* r=u*s */
#define _suNfc_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<10; ++_i){\
         _complex_mul((r).c[_i],(u).c[_k],(s).c[0]); ++_k;\
         int _j=1; for (; _j<8; ){ \
            _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
         } \
         _complex_mul_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
      }\
   } while(0) 

/* SU(N) matrix u^dagger times SU(N) vector s */
/* r=(u^dagger)*s */
#define _suNfc_inverse_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<10; ++_i){\
         _complex_mul_star((r).c[_i],(s).c[0],(u).c[_k]);\
         int _j=1; for (; _j<8; ){ \
            _k+=10; _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]); ++_j;\
            _k+=10; _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]); ++_j;\
            _k+=10; _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]); ++_j;\
            _k+=10; _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]); ++_j;\
         } \
         _k+=10; _complex_mul_star_assign((r).c[_i],(s).c[_j],(u).c[_k]); ++_j;\
         _k-=89;\
      }\
   } while(0) 

/* u=0 */
#define _suNfc_zero(u) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/* SU(N) matrix u times SU(N) vector s */
/* r=u*s */
#define _suNf_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<10; ++_i){\
         _complex_mulr((r).c[_i],(u).c[_k],(s).c[0]); ++_k;\
         int _j=1; for (; _j<8; ){ \
            _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
            _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
         } \
         _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_k; ++_j;\
      }\
   } while(0) 

/* SU(N) matrix u^dagger times SU(N) vector s */
/* r=(u^dagger)*s */
#define _suNf_inverse_multiply(r,u,s) \
   do { \
      int _i,_k=0;for (_i=0; _i<10; ++_i){\
         _complex_mulr((r).c[_i],(u).c[_k],(s).c[0]);\
         int _j=1; for (; _j<8; ){ \
            _k+=10; _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_j;\
            _k+=10; _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_j;\
            _k+=10; _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_j;\
            _k+=10; _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_j;\
         } \
         _k+=10; _complex_mulr_assign((r).c[_i],(u).c[_k],(s).c[_j]); ++_j;\
         _k-=89;\
      }\
   } while(0) 

/* u=0 */
#define _suNf_zero(u) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
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
#define _suNfc_dagger(u,v) \
   do { \
      int _i,_j,_n=0,_k=0;\
      for (_i=0; _i<10; ++_i){\
         _complex_star((u).c[_n],(v).c[_k]);\
         for (_j=0; _j<8; ){ \
            ++_n; _k+=10; _complex_star((u).c[_n],(v).c[_k]); ++_j;\
            ++_n; _k+=10; _complex_star((u).c[_n],(v).c[_k]); ++_j;\
            ++_n; _k+=10; _complex_star((u).c[_n],(v).c[_k]); ++_j;\
            ++_n; _k+=10; _complex_star((u).c[_n],(v).c[_k]); ++_j;\
         } \
         ++_n; _k+=10; _complex_star((u).c[_n],(v).c[_k]);\
         ++_n; _k-=89;\
      }\
   } while(0) 

/* u=v*w */
#define _suNfc_times_suNfc(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            _complex_mul((u).c[_n],(v).c[_k],(w).c[_l]);\
            for (_j=0; _j<8; ){ \
               ++_k; _l+=10; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; _l+=10; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; _l+=10; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; _l+=10; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
            } \
            ++_k; _l+=10; _complex_mul_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_n; _k-=9; _l-=89;\
         } _k+=10; _l=0;\
      }\
   } while(0) 

/* u=v*w^+ */
#define _suNfc_times_suNfc_dagger(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            _complex_mul_star((u).c[_n],(v).c[_k],(w).c[_l]);\
            for (_j=0; _j<8; ){ \
               ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
               ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]); ++_j;\
            } \
            ++_k; ++_l; _complex_mul_star_assign((u).c[_n],(v).c[_k],(w).c[_l]);\
            ++_n; _k-=9; ++_l;\
         } _k+=10; _l=0;\
      }\
   } while(0) 

/* u=v^+*w */
#define _suNfc_dagger_times_suNfc(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            _k=_y; _l=_i;\
            _complex_mul_star((u).c[_n],(w).c[_k],(v).c[_l]);\
            for (_j=0; _j<8; ){ \
               _k+=10; _l+=10; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]); ++_j;\
               _k+=10; _l+=10; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]); ++_j;\
               _k+=10; _l+=10; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]); ++_j;\
               _k+=10; _l+=10; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]); ++_j;\
            } \
            _k+=10; _l+=10; _complex_mul_star_assign((u).c[_n],(w).c[_k],(v).c[_l]);\
            ++_n;\
         }\
      }\
   } while(0) 

/* u=0 */
#define _suNfc_zero(u) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=1 */
#define _suNfc_unit(u) \
   do { \
      _suNfc_zero((u));\
      int _i,_n=0; for (_i=0; _i<8; ){\
         _complex_1((u).c[_n]); _n+=11; ++_i;\
         _complex_1((u).c[_n]); _n+=11; ++_i;\
         _complex_1((u).c[_n]); _n+=11; ++_i;\
         _complex_1((u).c[_n]); _n+=11; ++_i;\
      }\
      _complex_1((u).c[_n]); _n+=11;\
      _complex_1((u).c[_n]); _n+=11;\
   } while(0) 

/* u=-v */
#define _suNfc_minus(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
         _complex_minus((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=v */
#define _suNfc_copy(u,v) \
   (u)=(v)

/* u=r*v (u,v mat, r real) */
#define _suNfc_mul(u,r,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mulr((u).c[_i],(r),(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=r*v (u,v mat, r complex) */
#define _suNfc_mulc(u,r,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
         _complex_mul((u).c[_i],(r),(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u+=v */
#define _suNfc_add_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_add_assign((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* u-=v */
#define _suNfc_sub_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
         _complex_sub_assign((u).c[_i],(v).c[_i]); ++_i;\
      }\
   } while(0) 

/* k=| u |2 ) */
#define _suNfc_sqnorm(k,u) \
   do { \
      int _i;\
      (k)=0.;\
      for (_i=0; _i<100; ){\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
         (k)+=_complex_prod_re((u).c[_i],(u).c[_i]); ++_i;\
      }\
   } while(0) 

/* k=| 1 - u |2 ) */
#define _suNfc_sqnorm_m1(k,u) \
   do { \
      (k)=0.;\
      int _i,_j,_n=0,_l=0,_s2=0;\
      for(_i=0;_i<10;){\
         (k)+=_complex_prod_m1_re((u).c[_n],(u).c[_n]);\
         ++_n; _l+=10; \
         for(_j=_i+1;_j<10;++_j){\
            (k)+=_complex_prod_re((u).c[_n],(u).c[_n]);\
            (k)+=_complex_prod_re((u).c[_l],(u).c[_l]);\
            ++_n; _l+=10; \
         }\
         ++_i; _s2+=10; _n+=_i; _l-=99-_s2; \
      }\
   } while(0) 

/* k=Re Tr (u) */
#define _suNfc_trace_re(k,u) \
   do { \
      int _i,_n=0;\
      (k)=0.;\
      for (_i=0; _i<8; _i+=4){\
         (k)+=_complex_re((u).c[_n])+ \
              _complex_re((u).c[_n+11])+ \
              _complex_re((u).c[_n+22])+ \
              _complex_re((u).c[_n+33]);\
         _n+=44;\
      }\
      (k)+=_complex_re((u).c[_n])+ \
           _complex_re((u).c[_n+11]);\
   } while(0) 

/* k=Im Tr (u) */
#define _suNfc_trace_im(k,u) \
   do { \
      int _i,_n=0;\
      (k)=0.;\
      for (_i=0; _i<8; _i+=4){\
         (k)+=_complex_im((u).c[_n])+ \
              _complex_im((u).c[_n+11])+ \
              _complex_im((u).c[_n+22])+ \
              _complex_im((u).c[_n+33]);\
         _n+=44;\
      }\
      (k)+=_complex_im((u).c[_n])+ \
           _complex_im((u).c[_n+11]);\
   } while(0) 

/* This fuction computes the hmc force matrix */
/* this fuction accumulates on the original matrix u */
#define _suNfc_FMAT(u,s) \
   do { \
      int _i,_j,_n=0;\
      for (_i=0; _i<10; ++_i){\
         for (_j=0; _j<8; ){\
            _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
         }\
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
         ++_n; ++_j; \
         _complex_mul_star_assign((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
         _complex_mul_star_assign((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
         ++_n; ++_j; \
      }\
   } while(0) 

/* This fuction compute the hmc force matrix */
/* this fuction accumulates on the original matrix u */
#define _suNf_FMAT(u,s) \
   do { \
      int _i,_j,_n=0;\
      for (_i=0; _i<10; ++_i){\
         for (_j=0; _j<8; ){\
            _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
            _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
            _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
            ++_n; ++_j; \
         }\
         _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
         _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
         ++_n; ++_j; \
         _complex_mul_star_assign_re((u).c[_n],(s).c[0].c[_i],(s).c[2].c[_j]); \
         _complex_mul_star_assign_re((u).c[_n],(s).c[1].c[_i],(s).c[3].c[_j]); \
         ++_n; ++_j; \
      }\
   } while(0) 

/* u=1 */
#define _suNf_unit(u) \
   do { \
      _suNf_zero((u));\
      int _i,_n=0; for (_i=0; _i<8; ){\
         (u).c[_n]=1.; _n+=11; ++_i;\
         (u).c[_n]=1.; _n+=11; ++_i;\
         (u).c[_n]=1.; _n+=11; ++_i;\
         (u).c[_n]=1.; _n+=11; ++_i;\
      }\
      (u).c[_n]=1.; _n+=11;\
      (u).c[_n]=1.; _n+=11;\
   } while(0) 

/* u=v^dagger */
#define _suNf_dagger(u,v) \
   {\
      int _i,_j,_n=0,_k=0;\
      for (_i=0; _i<10; ++_i){\
         (u).c[_n]=(v).c[_k];\
         for (_j=0; _j<8; ){ \
            ++_n; _k+=10; (u).c[_n]=(v).c[_k]; ++_j;\
            ++_n; _k+=10; (u).c[_n]=(v).c[_k]; ++_j;\
            ++_n; _k+=10; (u).c[_n]=(v).c[_k]; ++_j;\
            ++_n; _k+=10; (u).c[_n]=(v).c[_k]; ++_j;\
         } \
         ++_n; _k+=10; (u).c[_n]=(v).c[_k];\
         ++_n; _k-=89;\
      }\
   }((void)0) 

/* u=v*w */
#define _suNf_times_suNf(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            (u).c[_n]=(v).c[_k]*(w).c[_l];\
            for (_j=0; _j<8; ){ \
               ++_k; _l+=10; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; _l+=10; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; _l+=10; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; _l+=10; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
            } \
            ++_k; _l+=10; (u).c[_n]+=(v).c[_k]*(w).c[_l];\
            ++_n; _k-=9; _l-=89;\
         } _k+=10; _l=0;\
      }\
   } while(0) 

/* u=v*w^t */
#define _suNf_times_suNf_dagger(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            (u).c[_n]=(v).c[_k]*(w).c[_l];\
            for (_j=0; _j<8; ){ \
               ++_k; ++_l; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; ++_l; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; ++_l; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
               ++_k; ++_l; (u).c[_n]+=(v).c[_k]*(w).c[_l]; ++_j;\
            } \
            ++_k; ++_l; (u).c[_n]+=(v).c[_k]*(w).c[_l];\
            ++_n; _k-=9; ++_l;\
         } _k+=10; _l=0;\
      }\
   } while(0) 

/* u=v^+*w */
#define _suNf_dagger_times_suNf(u,v,w) \
   do { \
      int _i,_y,_j,_n=0,_k=0,_l=0;\
      for (_i=0; _i<10; ++_i){\
         for (_y=0; _y<10; ++_y){\
            _k=_y; _l=_i;\
            (u).c[_n]=(w).c[_k]*(v).c[_l];\
            for (_j=0; _j<8; ){ \
               _k+=10; _l+=10; (u).c[_n]+=(w).c[_k]*(v).c[_l]; ++_j;\
               _k+=10; _l+=10; (u).c[_n]+=(w).c[_k]*(v).c[_l]; ++_j;\
               _k+=10; _l+=10; (u).c[_n]+=(w).c[_k]*(v).c[_l]; ++_j;\
               _k+=10; _l+=10; (u).c[_n]+=(w).c[_k]*(v).c[_l]; ++_j;\
            } \
            _k+=10; _l+=10; (u).c[_n]+=(w).c[_k]*(v).c[_l];\
            ++_n;\
         }\
      }\
   } while(0) 

/* u+=v */
#define _suNf_add_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]+=(v).c[_i]; ++_i;\
         (u).c[_i]+=(v).c[_i]; ++_i;\
         (u).c[_i]+=(v).c[_i]; ++_i;\
         (u).c[_i]+=(v).c[_i]; ++_i;\
      }\
   } while(0) 

/* u-=v */
#define _suNf_sub_assign(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]-=(v).c[_i]; ++_i;\
         (u).c[_i]-=(v).c[_i]; ++_i;\
         (u).c[_i]-=(v).c[_i]; ++_i;\
         (u).c[_i]-=(v).c[_i]; ++_i;\
      }\
   } while(0) 

/* u=r*v (u,v mat, r real) */
#define _suNf_mul(u,r,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]=(r)*(v).c[_i]; ++_i;\
         (u).c[_i]=(r)*(v).c[_i]; ++_i;\
         (u).c[_i]=(r)*(v).c[_i]; ++_i;\
         (u).c[_i]=(r)*(v).c[_i]; ++_i;\
      }\
   } while(0) 

/* u=v */
#define _suNf_copy(u,v) \
   (u)=(v)

/* k=Re Tr (u) */
#define _suNf_trace_re(k,u) \
   do { \
      int _i,_n=0;\
      (k)=0.;\
      for (_i=0; _i<8; _i+=4){\
         (k)+=(u).c[_n]+ \
              (u).c[_n+11]+ \
              (u).c[_n+22]+ \
              (u).c[_n+33];\
         _n+=44;\
      }\
      (k)+=(u).c[_n]+ \
           (u).c[_n+11];\
   } while(0) 

/* k= | u |2 ) */
#define _suNf_sqnorm(k,u) \
   do { \
      int _i;\
      (k)=0.;\
      for (_i=0; _i<100; ){\
         (k)+=(u).c[_i]*(u).c[_i]; ++_i;\
         (k)+=(u).c[_i]*(u).c[_i]; ++_i;\
         (k)+=(u).c[_i]*(u).c[_i]; ++_i;\
         (k)+=(u).c[_i]*(u).c[_i]; ++_i;\
      }\
   } while(0) 

/* u=-v */
#define _suNf_minus(u,v) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]=-(v).c[_i]; ++_i;\
         (u).c[_i]=-(v).c[_i]; ++_i;\
         (u).c[_i]=-(v).c[_i]; ++_i;\
         (u).c[_i]=-(v).c[_i]; ++_i;\
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
#define _suNfc_FMAT_zero(u) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
         _complex_0((u).c[_i]); ++_i;\
      }\
   } while(0) 

/* u=0 */
#define _suNf_FMAT_zero(u) \
   do { \
      int _i;for (_i=0; _i<100; ){\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
         (u).c[_i]=0.; ++_i;\
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
/* (output) v = suNfc_vector ; (input) in = suNfc_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNf_read_spinor_flt_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*10)*(stride); \
      (v).c[0]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[8]=((complex_flt*)(in))[__iz]; __iz+=(stride); \
      (v).c[9]=((complex_flt*)(in))[__iz]; \
   } while (0) 

#define _suNf_read_spinor_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*10)*(stride); \
      (v).c[0]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[1]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[2]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[3]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[4]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[5]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[6]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[7]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[8]=((complex*)(in))[__iz]; __iz+=(stride); \
      (v).c[9]=((complex*)(in))[__iz]; \
   } while (0) 

/* Write spinor field component to GPU memory */
/* (input) v = suNfc_vector ; (output) out = suNfc_spinor* */
/* (input) iy = site ; (input) x = 0..3 spinor component; */
#define _suNf_write_spinor_flt_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*10)*(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[7]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[8]; __iz+=(stride); \
      ((complex_flt*)(out))[__iz]=(v).c[9]; \
   } while (0) 

#define _suNf_write_spinor_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*10)*(stride); \
      ((complex*)(out))[__iz]=(v).c[0]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[1]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[2]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[3]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[4]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[5]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[6]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[7]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[8]; __iz+=(stride); \
      ((complex*)(out))[__iz]=(v).c[9]; \
   } while (0) 

/* Read an suN algebra vector from GPU memory */
/* (output) v = suN_algebra_vector ; (input) in = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_av_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      (v).c[34]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[35]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[36]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[37]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[38]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[39]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[40]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[41]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[42]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[43]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[44]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[45]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[46]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[47]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[48]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[49]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[50]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[51]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[52]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[53]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[54]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[55]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[56]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[57]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[58]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[59]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[60]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[61]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[62]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[63]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[64]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[65]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[66]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[67]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[68]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[69]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[70]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[71]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[72]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[73]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[74]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[75]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[76]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[77]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[78]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[79]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[80]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[81]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[82]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[83]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[84]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[85]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[86]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[87]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[88]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[89]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[90]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[91]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[92]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[93]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[94]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[95]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[96]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[97]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[98]=((float*)(in))[__iz]; \
   } while (0) 

#define _suNf_av_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      (v).c[34]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[35]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[36]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[37]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[38]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[39]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[40]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[41]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[42]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[43]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[44]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[45]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[46]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[47]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[48]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[49]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[50]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[51]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[52]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[53]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[54]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[55]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[56]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[57]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[58]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[59]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[60]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[61]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[62]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[63]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[64]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[65]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[66]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[67]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[68]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[69]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[70]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[71]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[72]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[73]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[74]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[75]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[76]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[77]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[78]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[79]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[80]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[81]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[82]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[83]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[84]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[85]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[86]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[87]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[88]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[89]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[90]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[91]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[92]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[93]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[94]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[95]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[96]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[97]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[98]=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN algebra vector to GPU memory */
/* (input) v = suN_algebra_vector ; (output) out = suN_algebra_vector* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_av_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      ((float*)(out))[__iz]=(v).c[34]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[35]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[36]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[37]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[38]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[39]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[40]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[41]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[42]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[43]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[44]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[45]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[46]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[47]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[48]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[49]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[50]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[51]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[52]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[53]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[54]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[55]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[56]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[57]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[58]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[59]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[60]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[61]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[62]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[63]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[64]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[65]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[66]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[67]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[68]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[69]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[70]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[71]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[72]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[73]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[74]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[75]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[76]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[77]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[78]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[79]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[80]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[81]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[82]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[83]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[84]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[85]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[86]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[87]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[88]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[89]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[90]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[91]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[92]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[93]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[94]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[95]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[96]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[97]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[98]; \
   } while (0) 

#define _suNf_av_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      ((double*)(out))[__iz]=(v).c[34]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[35]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[36]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[37]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[38]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[39]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[40]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[41]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[42]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[43]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[44]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[45]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[46]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[47]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[48]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[49]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[50]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[51]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[52]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[53]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[54]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[55]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[56]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[57]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[58]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[59]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[60]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[61]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[62]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[63]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[64]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[65]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[66]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[67]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[68]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[69]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[70]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[71]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[72]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[73]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[74]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[75]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[76]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[77]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[78]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[79]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[80]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[81]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[82]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[83]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[84]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[85]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[86]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[87]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[88]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[89]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[90]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[91]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[92]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[93]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[94]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[95]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[96]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[97]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[98]; \
   } while (0) 

/* Mul_add_assign on a suN algebra vector on GPU  */
/* (in/out) v = suN_algebra_vector* ; (input) in = suN_algebra_vector */
/* (input) iy = site ; (input) x = 0..3 direction; (input) r = real */
#define _algebra_vector_mul_add_assign_gpu_f_flt(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      ((float*)(v))[__iz]+=(in).c[34]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[35]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[36]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[37]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[38]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[39]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[40]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[41]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[42]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[43]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[44]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[45]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[46]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[47]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[48]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[49]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[50]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[51]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[52]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[53]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[54]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[55]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[56]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[57]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[58]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[59]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[60]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[61]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[62]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[63]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[64]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[65]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[66]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[67]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[68]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[69]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[70]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[71]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[72]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[73]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[74]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[75]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[76]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[77]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[78]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[79]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[80]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[81]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[82]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[83]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[84]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[85]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[86]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[87]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[88]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[89]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[90]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[91]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[92]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[93]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[94]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[95]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[96]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[97]*(r); __iz+=(stride); \
      ((float*)(v))[__iz]+=(in).c[98]*(r); \
   } while (0) 

#define _algebra_vector_mul_add_assign_gpu_f(stride,v,iy,x,r,in) \
   do {  \
      int __iz=(iy)+((x)*99)*(stride); \
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
      ((double*)(v))[__iz]+=(in).c[34]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[35]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[36]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[37]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[38]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[39]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[40]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[41]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[42]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[43]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[44]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[45]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[46]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[47]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[48]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[49]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[50]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[51]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[52]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[53]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[54]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[55]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[56]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[57]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[58]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[59]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[60]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[61]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[62]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[63]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[64]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[65]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[66]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[67]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[68]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[69]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[70]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[71]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[72]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[73]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[74]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[75]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[76]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[77]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[78]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[79]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[80]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[81]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[82]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[83]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[84]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[85]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[86]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[87]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[88]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[89]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[90]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[91]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[92]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[93]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[94]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[95]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[96]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[97]*(r); __iz+=(stride); \
      ((double*)(v))[__iz]+=(in).c[98]*(r); \
   } while (0) 

/* Read an suN matrix from GPU memory */
/* (output) v = suN ; (input) in = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*100)*(stride); \
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
      (v).c[34]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[35]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[36]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[37]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[38]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[39]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[40]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[41]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[42]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[43]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[44]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[45]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[46]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[47]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[48]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[49]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[50]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[51]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[52]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[53]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[54]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[55]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[56]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[57]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[58]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[59]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[60]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[61]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[62]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[63]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[64]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[65]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[66]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[67]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[68]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[69]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[70]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[71]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[72]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[73]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[74]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[75]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[76]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[77]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[78]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[79]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[80]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[81]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[82]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[83]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[84]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[85]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[86]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[87]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[88]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[89]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[90]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[91]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[92]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[93]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[94]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[95]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[96]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[97]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[98]=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[99]=((float*)(in))[__iz]; \
   } while (0) 

#define _suNf_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*100)*(stride); \
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
      (v).c[34]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[35]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[36]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[37]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[38]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[39]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[40]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[41]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[42]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[43]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[44]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[45]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[46]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[47]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[48]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[49]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[50]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[51]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[52]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[53]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[54]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[55]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[56]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[57]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[58]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[59]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[60]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[61]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[62]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[63]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[64]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[65]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[66]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[67]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[68]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[69]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[70]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[71]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[72]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[73]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[74]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[75]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[76]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[77]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[78]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[79]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[80]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[81]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[82]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[83]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[84]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[85]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[86]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[87]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[88]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[89]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[90]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[91]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[92]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[93]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[94]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[95]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[96]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[97]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[98]=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[99]=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN matrix to GPU memory */
/* (input) v = suN ; (output) out = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNf_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*100)*(stride); \
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
      ((float*)(out))[__iz]=(v).c[34]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[35]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[36]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[37]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[38]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[39]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[40]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[41]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[42]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[43]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[44]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[45]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[46]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[47]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[48]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[49]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[50]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[51]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[52]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[53]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[54]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[55]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[56]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[57]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[58]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[59]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[60]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[61]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[62]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[63]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[64]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[65]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[66]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[67]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[68]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[69]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[70]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[71]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[72]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[73]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[74]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[75]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[76]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[77]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[78]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[79]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[80]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[81]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[82]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[83]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[84]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[85]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[86]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[87]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[88]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[89]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[90]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[91]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[92]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[93]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[94]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[95]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[96]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[97]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[98]; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[99]; \
   } while (0) 

#define _suNf_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*100)*(stride); \
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
      ((double*)(out))[__iz]=(v).c[34]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[35]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[36]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[37]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[38]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[39]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[40]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[41]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[42]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[43]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[44]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[45]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[46]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[47]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[48]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[49]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[50]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[51]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[52]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[53]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[54]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[55]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[56]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[57]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[58]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[59]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[60]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[61]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[62]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[63]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[64]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[65]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[66]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[67]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[68]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[69]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[70]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[71]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[72]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[73]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[74]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[75]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[76]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[77]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[78]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[79]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[80]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[81]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[82]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[83]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[84]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[85]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[86]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[87]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[88]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[89]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[90]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[91]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[92]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[93]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[94]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[95]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[96]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[97]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[98]; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[99]; \
   } while (0) 

/* Read an suN matrix from GPU memory */
/* (output) v = suN ; (input) in = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNfc_flt_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*200)*(stride); \
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
      (v).c[35].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[36].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[36].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[37].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[37].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[38].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[38].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[39].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[39].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[40].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[40].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[41].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[41].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[42].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[42].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[43].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[43].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[44].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[44].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[45].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[45].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[46].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[46].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[47].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[47].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[48].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[48].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[49].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[49].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[50].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[50].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[51].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[51].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[52].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[52].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[53].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[53].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[54].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[54].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[55].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[55].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[56].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[56].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[57].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[57].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[58].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[58].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[59].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[59].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[60].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[60].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[61].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[61].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[62].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[62].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[63].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[63].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[64].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[64].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[65].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[65].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[66].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[66].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[67].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[67].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[68].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[68].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[69].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[69].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[70].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[70].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[71].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[71].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[72].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[72].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[73].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[73].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[74].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[74].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[75].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[75].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[76].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[76].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[77].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[77].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[78].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[78].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[79].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[79].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[80].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[80].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[81].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[81].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[82].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[82].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[83].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[83].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[84].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[84].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[85].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[85].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[86].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[86].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[87].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[87].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[88].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[88].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[89].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[89].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[90].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[90].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[91].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[91].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[92].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[92].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[93].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[93].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[94].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[94].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[95].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[95].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[96].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[96].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[97].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[97].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[98].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[98].im=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[99].re=((float*)(in))[__iz]; __iz+=(stride); \
      (v).c[99].im=((float*)(in))[__iz]; \
   } while (0) 

#define _suNfc_read_gpu(stride,v,in,iy,x) \
   do {  \
      int __iz=(iy)+((x)*200)*(stride); \
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
      (v).c[35].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[36].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[36].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[37].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[37].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[38].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[38].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[39].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[39].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[40].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[40].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[41].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[41].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[42].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[42].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[43].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[43].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[44].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[44].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[45].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[45].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[46].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[46].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[47].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[47].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[48].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[48].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[49].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[49].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[50].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[50].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[51].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[51].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[52].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[52].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[53].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[53].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[54].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[54].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[55].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[55].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[56].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[56].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[57].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[57].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[58].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[58].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[59].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[59].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[60].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[60].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[61].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[61].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[62].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[62].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[63].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[63].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[64].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[64].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[65].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[65].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[66].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[66].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[67].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[67].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[68].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[68].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[69].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[69].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[70].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[70].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[71].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[71].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[72].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[72].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[73].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[73].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[74].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[74].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[75].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[75].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[76].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[76].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[77].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[77].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[78].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[78].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[79].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[79].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[80].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[80].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[81].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[81].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[82].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[82].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[83].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[83].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[84].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[84].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[85].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[85].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[86].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[86].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[87].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[87].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[88].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[88].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[89].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[89].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[90].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[90].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[91].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[91].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[92].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[92].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[93].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[93].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[94].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[94].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[95].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[95].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[96].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[96].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[97].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[97].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[98].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[98].im=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[99].re=((double*)(in))[__iz]; __iz+=(stride); \
      (v).c[99].im=((double*)(in))[__iz]; \
   } while (0) 

/* Write an suN matrix to GPU memory */
/* (input) v = suN ; (output) out = suN* */
/* (input) iy = site ; (input) x = 0..3 direction; */
#define _suNfc_flt_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*200)*(stride); \
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
      ((float*)(out))[__iz]=(v).c[35].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[36].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[36].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[37].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[37].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[38].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[38].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[39].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[39].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[40].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[40].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[41].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[41].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[42].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[42].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[43].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[43].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[44].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[44].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[45].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[45].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[46].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[46].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[47].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[47].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[48].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[48].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[49].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[49].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[50].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[50].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[51].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[51].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[52].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[52].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[53].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[53].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[54].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[54].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[55].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[55].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[56].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[56].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[57].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[57].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[58].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[58].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[59].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[59].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[60].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[60].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[61].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[61].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[62].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[62].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[63].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[63].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[64].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[64].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[65].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[65].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[66].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[66].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[67].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[67].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[68].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[68].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[69].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[69].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[70].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[70].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[71].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[71].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[72].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[72].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[73].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[73].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[74].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[74].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[75].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[75].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[76].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[76].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[77].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[77].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[78].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[78].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[79].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[79].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[80].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[80].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[81].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[81].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[82].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[82].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[83].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[83].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[84].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[84].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[85].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[85].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[86].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[86].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[87].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[87].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[88].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[88].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[89].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[89].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[90].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[90].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[91].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[91].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[92].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[92].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[93].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[93].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[94].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[94].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[95].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[95].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[96].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[96].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[97].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[97].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[98].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[98].im; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[99].re; __iz+=(stride); \
      ((float*)(out))[__iz]=(v).c[99].im; \
   } while (0) 

#define _suNfc_write_gpu(stride,v,out,iy,x) \
   do {  \
      int __iz=(iy)+((x)*200)*(stride); \
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
      ((double*)(out))[__iz]=(v).c[35].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[36].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[36].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[37].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[37].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[38].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[38].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[39].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[39].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[40].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[40].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[41].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[41].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[42].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[42].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[43].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[43].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[44].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[44].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[45].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[45].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[46].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[46].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[47].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[47].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[48].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[48].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[49].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[49].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[50].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[50].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[51].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[51].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[52].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[52].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[53].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[53].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[54].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[54].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[55].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[55].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[56].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[56].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[57].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[57].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[58].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[58].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[59].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[59].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[60].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[60].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[61].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[61].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[62].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[62].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[63].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[63].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[64].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[64].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[65].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[65].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[66].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[66].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[67].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[67].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[68].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[68].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[69].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[69].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[70].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[70].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[71].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[71].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[72].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[72].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[73].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[73].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[74].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[74].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[75].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[75].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[76].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[76].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[77].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[77].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[78].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[78].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[79].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[79].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[80].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[80].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[81].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[81].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[82].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[82].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[83].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[83].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[84].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[84].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[85].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[85].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[86].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[86].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[87].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[87].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[88].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[88].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[89].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[89].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[90].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[90].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[91].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[91].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[92].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[92].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[93].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[93].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[94].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[94].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[95].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[95].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[96].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[96].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[97].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[97].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[98].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[98].im; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[99].re; __iz+=(stride); \
      ((double*)(out))[__iz]=(v).c[99].im; \
   } while (0) 


#endif
