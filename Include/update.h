/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"
#include "inverters.h"
#include "rational_functions.h"
#include "glueballs.h"

void test_staples();


void update(double *beta, int nhb, int nor);
void random_su2(double rho, double s[]);


/* functions and structures for the MRE algorithm */
typedef struct
{
	spinor_field *s[2];
	int num[2];
	int max;
	int init;
} mre_par;


typedef struct
{
	int id;
	int n_pf;
	spinor_field *pf;
	double mass;
	rational_app *ratio;
	double inv_err2;
	suNg_av_field **momenta;
} force_rhmc_par;

typedef struct
{
	int id;
	int n_pf;
	spinor_field *pf;
	int hasenbusch;
	double mass;
	double b;
	double mu;
	double inv_err2, inv_err2_flt;
	mre_par mpar;
	int logdet;
	suNg_av_field **momenta;
} force_hmc_par;

typedef struct
{
	double beta;
	double c0;
	double c1;
	suNg_av_field **momenta;
} force_gauge_par;

typedef struct
{
	double mass;
	double lambda;
	suNg_scalar_field **momenta;
	suNg_av_field **g_momenta;
} force_scalar_par;

typedef struct
{
	double gamma;
} force_auxfield_par;

typedef struct
{
	suNg_field **field;
	suNg_av_field **momenta;
} field_gauge_par;

typedef struct
{
	suNg_scalar_field **field;
	suNg_scalar_field **momenta;
} field_scalar_par;

void update_gauge_field(double, void *);
void update_auxfields(double, void *);

void update_scalar_field(double, void *);
void force_scalar(double, void *);


void fermion_force_begin();
void fermion_force_end(double dt, suNg_av_field *);
void force_fermion_core(spinor_field *, spinor_field *, int, double, double);
void force_fermion_core_taylor(spinor_field *, spinor_field *, int, double, double);
void force_clover_logdet(double, double);
void force_clover_fermion(spinor_field *Xs, spinor_field *Ys, double residue);
void force_clover_fermion_taylor(spinor_field *Xs, spinor_field *Ys, double residue);

void force_hmc(double, void *);
void force_hmc_tm(double, void *);
void force_rhmc(double, void *);
void force0(double, void *);
void force_hmc_auxfields(double, void *); //Force from a four_fermion monomial
void force_hmc_ff(double, void *);		  //Force from a HMC_ff or Hasenbusch_ff monomial


/* For the fermion force ? */
void corret_pf_dist_hmc();
void calc_one_force(int n_force);

#include "monomials.h"

typedef struct _integrator_par
{
	int nsteps;
	int nmon;
	const monomial **mon_list;
	void (*integrator)(double, struct _integrator_par *);
	struct _integrator_par *next;
	int level;
} integrator_par;

void leapfrog_multistep(double tlen, integrator_par *int_par);
void O2MN_multistep(double tlen, integrator_par *int_par);
void O4MN_multistep(double tlen, integrator_par *int_par);

typedef struct _ghmc_par
{

	/* integrator */
	integrator_par *integrator;
	double tlen;
	double csw;
	double rho_s;
	double rho_t;

	/* Fermion Theta angles */
	double theta[4];

	/* SF stuff */
	double SF_zf;
	double SF_ds;
	int SF_sign;
	double SF_ct;
	int SF_background;

} ghmc_par;

void free_ghmc();

/* stout smearing */

/* local action */
typedef enum
{
	NEW = 1,
	DELTA = 2
} local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */



/* Utility functions for four fermion interactions */

#endif
