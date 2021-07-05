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
#include "inverters.h"
#include "geometry.h"
#include <stdlib.h>

/*SUN exp matrix*/

void ExpX(double dt, suNg_algebra_vector *h, suNg *u);


typedef struct
{
  double gauge_boundary_improvement_cs;
  double gauge_boundary_improvement_ct;
  double chiSF_boundary_improvement_ds;
  double fermion_twisting_theta[4];
  int SF_BCs;
  suNg gauge_boundary_up;
  suNg gauge_boundary_dn;
} BCs_pars_t;

void init_BCs(BCs_pars_t *pars);
void init_plaq_open_BCs(double *plaq_weight, double *rect_weight, double ct, double cs);

void free_BCs();
void apply_BCs_on_represented_gauge_field();
void apply_BCs_on_fundamental_gauge_field();
void apply_BCs_on_spinor_field(spinor_field *sp);

#if defined(GAUGE_SPN) && defined(REPR_FUNDAMENTAL)
void apply_BCs_on_clover_term(suNffull_field *);
#else
void apply_BCs_on_clover_term(suNfc_field *);
#endif


void SF_classical_solution();


/*Global shift for fields, the routine accepts also NULL entries in which case it does nothing*/

void cross_prod(suNg_vector *v1, suNg_vector *v2, suNg_vector *v3);
void cross_prod_flt(suNg_vector_flt *v1, suNg_vector_flt *v2, suNg_vector_flt *v3);

#ifndef GAUGE_SON
#else
int project_to_suNg_real(suNg *out, suNg *in);
void diag_hmat(suNg *hmat, double *dag);
#endif

void assign_u2ud(void);
void assign_ud2u(void);
void assign_ud2u_f(void);



/* use power method to find max eigvalue of H2 */

/* EVA preconditioning */
typedef struct _eva_prec
{
  /* EVA parameters */
  int nevt;      /* search space dimension */
  int nev;       /* number of accurate eigenvalues */
  int kmax;      /* max degree of polynomial */
  int maxiter;   /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
} eva_prec;

/* HYP smearing */
void spatialHYP_smearing(suNg_field *out, suNg_field *in, double weight[3]);

/* Timing */
#include <sys/time.h>
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

/* CINFO */
void print_compiling_info();
void print_compiling_info_short();

/* Spatial Trasformations*/
#ifdef MAIN_PROGRAM
int *active_slices_list = NULL;
int *glbT_to_active_slices = NULL;
int n_active_slices;
#else
extern int *active_slices_list;
extern int *glbT_to_active_slices;
extern int n_active_slices;
#endif

/* Spatial blocking */
typedef enum
{
  NEW_SBLK = 1,
  CONT_SBLK = 0
} eval_spat_block;


/* Spatial rotation*/

/* Spatial APE smearing*/

/* Workspace database*/
void free_wrk_space();

#endif
