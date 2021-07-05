/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "hr_complex.h"
#include "spinor_field.h"


typedef void (*spinor_operator)(spinor_field *out, spinor_field *in);
typedef void (*spinor_operator_flt)(spinor_field_flt *out, spinor_field_flt *in);
typedef void (*spinor_operator_m)(spinor_field *out, spinor_field *in, double m);

typedef struct _mshift_par {
   int n; /* number of shifts */
   double *shift;
   double err2; /* relative error of the solutions */
   int max_iter; /* maximum number of iterations: 0 => infinity */
	 void *add_par; /* additional parameters for specific inverters */
} mshift_par;



// We might want to add a "trial" spinor to the argument list
typedef int (*inverter_ptr)(mshift_par* par, spinor_operator M, spinor_field *in, spinor_field *out);


/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int cg_mshift_flt(mshift_par *par, spinor_operator M, spinor_operator_flt F, spinor_field *in, spinor_field *out);


typedef struct {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
  double err2_flt; /* maximum error on the solutions */
  int max_iter_flt; /* maximum number of iterations: 0 => infinity */
} g5QMR_fltacc_par;

int g5QMR_mshift_trunc(mshift_par *par, int trunc_iter, spinor_operator M, spinor_field *in, spinor_field *out_trunc, spinor_field *out);


typedef struct _MINRES_par {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_par;


/*
int BiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
int HBiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
*/









#endif
