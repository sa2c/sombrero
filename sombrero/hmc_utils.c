/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "logger.h"
#include "random.h"
#include "representation.h"
#include "update.h"
#include "memory.h"
#include "utils.h"
#include "communications.h"
#include "clover_tools.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif



/* Initialize the Monte Carlo.
 * This performs the following operations:
 * 1) set the starting gauge field
 * 2) init the hmc update
 */
int init_mc() {

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  u_gauge_f_flt=alloc_gfield_f_flt(&glattice);

#ifdef WITH_CLOVER
	clover_init(1.0);
#endif

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 0
  };
  init_BCs(&BCs_pars);

  /* init gauge field */
  unit_u(u_gauge);
  
  apply_BCs_on_fundamental_gauge_field(); 
  represent_gauge_field();

  return 0;
}
    
/* clean up memory */
int end_mc() {
  free_BCs();
  
  /* free memory */
  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  return 0;
}
