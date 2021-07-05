/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *
 * All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include "monomials.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"
#include "observables.h"

/* State quantities for HMC */
static suNg_field *u_gauge_old = NULL;
static suNg_scalar_field *u_scalar_old = NULL;
static scalar_field *ff_sigma_old = NULL;
static scalar_field *ff_pi_old = NULL;

static scalar_field *la = NULL; /* local action field for Metropolis test */

static ghmc_par update_par;
static int init = 0;

  //#ifdef ROTATED_SF
  //#endif

void free_ghmc()
{
  if (!init)
    return;

  /* free momenta */
  if (u_gauge_old != NULL)
  {
    free_gfield(u_gauge_old);
    u_gauge_old = NULL;
  }
  if (u_scalar_old != NULL)
  {
    free_scalar_field(u_scalar_old);
    u_scalar_old = NULL;
  }
  if (suN_momenta != NULL)
  {
    free_avfield(suN_momenta);
    suN_momenta = NULL;
  }
  if (scalar_momenta != NULL)
  {
    free_scalar_field(scalar_momenta);
    scalar_momenta = NULL;
  }
  if (la != NULL)
  {
    free_sfield(la);
    la = NULL;
  }

  /*Free integrator */
  integrator_par *ip = update_par.integrator;
  while (ip != NULL)
  {
    update_par.integrator = ip->next;
    free(ip->mon_list);
    free(ip);
    ip = update_par.integrator;
  }
  update_par.integrator = NULL;

  //free_force_hmc();
  init = 0;
  lprintf("HMC", 0, "Memory deallocated.\n");
}



/*
### this is essentially a copy of the function update_ghmc which flips momenta.
### this is to allow for a reversibility test.
### this might have to be changed  if update_ghmc is modified.
*/


#ifdef MEASURE_FORCEHMC
/*Functions to check forces */
void corret_pf_dist_hmc()
{
  /* init monomials */
  for (int i = 0; i < num_mon(); ++i)
  {
    const monomial *m = mon_n(i);
    m->init_traj(m);
  }

  /* generate new momenta */
  lprintf("HMC", 30, "Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(suN_momenta);

  /* generate new pseudofermions */
  for (int i = 0; i < num_mon(); ++i)
  {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("HMC", 30, "Computing action density...\n");
  local_hmc_action(NEW, la, suN_momenta, scalar_momenta);

  /* correct pseudofermion distribution */
  for (int i = 0; i < num_mon(); ++i)
  {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }
}

void calc_one_force(int n_force)
{
  integrator_par *ip = update_par.integrator;
  for (;;)
  {
    error(ip == NULL, 1, "calc_one_force", "Error in force index\n");
    for (int n = 0; n < ip->nmon; n++)
    {
      const monomial *m = ip->mon_list[n];
      if (m->data.id == n_force)
      {
        m->force_f(1, m->force_par);
        return;
      }
    }
    ip = ip->next;
  }
}
#endif
