/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *
 * All rights reserved.                                                      *
 \***************************************************************************/

/*******************************************************************************
 *
 * File process_init.c
 *
 * Inizialization of geometry structures
 *
 *******************************************************************************/

#include "global.h"
#include "error.h"
#include "logger.h"
#include "hr_omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "io.h"
#include "random.h"
#include "setup.h"
#include "memory.h"
#include "representation.h"
#include "clover_tools.h"
#include "clover_exp.h"

/* setup_process
 * Assign a unique RID, PID to each process and setup
 * communications as necessary
 *
 * return codes:
 * 0 => success
 *
 * OUTPUT:
 * MPI: GLB_COMM, RID, PID, WORLD_SIZE
 */

static char input_filename[256] = "input_file";
static char output_filename[256] = "out_0";
static char error_filename[256] = "err_0";



static int setup_level = 0;


#ifndef REPR_FUNDAMENTAL
#endif
#ifndef ALLOCATE_REPR_GAUGE_FIELD
#endif
#ifdef DPHI_FLT
#endif
#ifdef WITH_SMEARING
#endif
#if defined(WITH_CLOVER)|| defined(WITH_EXPCLOVER)
#endif

#ifdef WITH_MPI
#else
#endif
#ifdef _OPENMP
#endif
#ifdef GAUGE_SON
#endif


/* this function is intended to clean up before process ending
 *
 * return codes:
 * 0 => success
 */
int finalize_process()
{

  free_ghmc();

  free_BCs();

  free_wrk_space();

  /* free memory */
  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif
  if (u_scalar != NULL)
    free_scalar_field(u_scalar);

  if (u_gauge_f_flt != NULL)
    free_gfield_f_flt(u_gauge_f_flt);

  free_geometry_mpi_eo();

  lprintf("SYSTEM", 0, "Process finalized.\n");

#ifdef WITH_MPI
  /* MPI variables */
  int init;
  MPI_Initialized(&init);
  if (init)
    MPI_Finalize();
#endif

  return 0;
}

/* setup_replicas
 * Split MPI_COMM_WORLD into replica communicators GLB_COMM
 *
 * return codes:
 * 0 => success
 *
 * AFFECTS THE GLOBAL VARIABLES: GLB_COMM, RID, PID, WORLD_SIZE
 */
#ifdef WITH_MPI
#endif //ifdef WITH_MPI
