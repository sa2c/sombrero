/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File utils.h
*
* Some useful functions
*
*******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "suN_types.h"
#include "spinor_field.h"
#include "geometry.h"

void ExpX(double dt, suNg_algebra_vector *h, suNg *u);

void vector_star(suNg_vector*,suNg_vector*);


typedef struct {
  double gauge_boundary_improvement_cs;
  double gauge_boundary_improvement_ct;
  double chiSF_boundary_improvement_ds;
  double fermion_twisting_theta[4];
  int SF_BCs;
  suNg gauge_boundary_up;
  suNg gauge_boundary_dn;
} BCs_pars_t;

void init_BCs(BCs_pars_t *pars);
void free_BCs();
void apply_BCs_on_represented_gauge_field();
void apply_BCs_on_fundamental_gauge_field();
void apply_BCs_on_momentum_field(suNg_av_field *force);
void apply_BCs_on_spinor_field(spinor_field *sp);
void apply_BCs_on_spinor_field_flt(spinor_field_flt *sp);
void apply_background_field_zdir(suNg_field* V,double Q,int n);

#if defined(GAUGE_SPN) && defined(REPR_FUNDAMENTAL)
void apply_BCs_on_clover_term(suNffull_field*);
#else
void apply_BCs_on_clover_term(suNfc_field*);
#endif

void project_to_suNg(suNg *u);
void project_to_suNg_flt(suNg_flt *u);
void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling);

#ifndef GAUGE_SON
void ludcmp(complex* a, int* indx, double* d,int N);
void lubksb(complex* a, int* indx, complex* b,int N);
void inv_suNg(suNg* a);
void det_suNg(complex* res, suNg* a);
#else
int project_to_suNg_real(suNg *out, suNg *in);
void det_suNg(double* res, suNg *a);
void diag_hmat(suNg *hmat, double *dag);
#endif

void assign_u2ud(void);
void assign_ud2u(void);
void assign_ud2u_f(void);

/* void assign_s2sd(int len, suNf_spinor *out, suNf_spinor_flt *in); */
/* void assign_sd2s(int len, suNf_spinor_flt *out, suNf_spinor *in); */

void assign_s2sd(spinor_field *out, spinor_field_flt *in);
void assign_sd2s(spinor_field_flt *out, spinor_field *in);

/* Timing */
#include <sys/time.h>
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);

#endif
