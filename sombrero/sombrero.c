/****************************************************************************
* Copyright (c) 2008, Claudio Pica                                          *
* Copyright (c) 2019, Jarno Rantaharju                                      *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
 *
 * Main HMC program
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dirac.h"
#include "global.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "logger.h"
#include "memory.h"
#include "random.h"
#include "sombrero.h"
#include "update.h"

#include "counters/communication_count.h"
#include "counters/flop_count.h"
#include "counters/memory_count.h"

#if defined __SSE3__ && defined __GNUC__
#include <pmmintrin.h>
#include <xmmintrin.h>
#endif

static int cg_test(spinor_field *in, spinor_field *out, int iterations);
int init_mc();
int end_mc();

/* Sort the directions based on size */
static void sort_dirs(int *sorted, int *vol) {
  sorted[0] = 0;
  sorted[1] = 1;
  sorted[2] = 2;
  sorted[3] = 3;
  for (int d = 0; d < 3; d++)
    for (int d2 = d + 1; d2 < 4; d2++) {
      if (vol[sorted[d2]] > vol[sorted[d]]) {
        int dt = sorted[d];
        sorted[d] = sorted[d2];
        sorted[d2] = dt;
      }
    }
}

/* Factorize an integer */
static void get_factors(int n, int *factors, int *nfactors) {
  *nfactors = 0;
  /* Check for factors and build a list */
  for (int f = 2; f <= n; f++) {
    while (n % f == 0) {
      factors[*nfactors] = f;
      n /= f;
      (*nfactors)++;
    }
  }
}

/* Multiply given local lattice size by a factor */
static void weak_scale_lattice(int x[4], int scale) {
  /* Get factors */
  int factors[scale];
  int n_factors;
  get_factors(scale, factors, &n_factors);

  /* Loop over factors in reverse order (largest first) and multiply the lattice
   * size */
  int dir = 0;
  for (int f = n_factors - 1; f >= 0; f--) {
    int factor = factors[f];
    x[dir] *= factor;
    dir = (dir + 1) % 4;
  }
}

/* Partition the lattice accross the MPI ranks */
static void split_lattice(int np[4]) {
  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  error(GLB_T * GLB_X * GLB_Y * GLB_Z % n_nodes, 1, "main [sombrero.c]",
        " Lattice volume must be divisible by the number of ranks\n");

  int d = 0, dirs[4] = {0, 1, 2, 3};
  int local[4] = {GLB_T, GLB_X, GLB_Y, GLB_Z};

  /* Get factors */
  int factors[n_nodes];
  int n_factors;
  get_factors(n_nodes, factors, &n_factors);

  /* Loop over factors (in reverse order, largest first) */
  sort_dirs(dirs, local);
  for (int f = n_factors - 1; f >= 0; f--) {
    int fac = factors[f];
    if (local[dirs[d]] % fac == 0) {
      np[dirs[d]] *= fac;
      local[dirs[d]] /= fac;
      n_nodes /= fac;
      sort_dirs(dirs, local);
      d = 0;
    } else {
      d = d + 1;
      f = f + 1;
      if (d > 3) {
        char message[200];
        sprintf(message, " Lattice volume cannot be divided across given "
                         "number of nodes.\n");
        sprintf(message,
                "%sAt least on factor of %d in the number of ranks cannot be "
                "used to divide any of the lattice directions.\n",
                message, fac);
        sprintf(message,
                "%sIn the default cases, use only factors of 2 and 3.\n",
                message);
        error(1, 1, "main [sombrero.c]", message);
      }
    }
  }
  for (int d = 0; d < 4; d++) {
    error(local[d] < 4, 1, "main [sombrero.c]",
          " Too small local lattice in direction.\nToo many MPI ranks or too "
          "small problem size.\n");
  }
}

static void read_cmdline_sombrero(int argc, char *argv[]) {
  int i, as = 0, aw = 0, ap = 0, av = 0, al = 0, requested = 1;

  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-s") == 0) {
      as = i + 1;
      requested += 2;
    } else if (strcmp(argv[i], "-w") == 0) {
      aw = i + 1;
      requested += 1;
    } else if (strcmp(argv[i], "-p") == 0) {
      ap = i + 1;
      requested += 2;
    } else if (strcmp(argv[i], "-l") == 0) {
      al = i + 1;
      requested += 2;
    } else if (strcmp(argv[i], "-v") == 0) {
      av = i + 1;
      requested += 2;
    }
  }

  error(argc < requested, 1, "read_cmdline_sombrero [sombrero.c]",
        "Arguments: [ -w ] [ -s small | medium | large | very_large ]");

  int loglevel = 10;
  if (av != 0) {
    if (strcmp(argv[av], "result") == 0) {
      loglevel = 0;
    } else if (strcmp(argv[av], "verbose") == 0) {
      loglevel = 50;
    } else {
      error(1, 1, "read_cmdline_sombrero [sombrero.c]",
            " Unrecognized verbosity level (use -v result or -v verbose)\n");
    }
  }

  /* set logger levels */
  input_logger logger_var;
  logger_var.def_log_lvl = loglevel;
  logger_var.inverter_log_lvl = loglevel;
  logger_var.forcestat_log_lvl = loglevel;

  /* run logger setup */
  logger_set_input(&logger_var);

  if (PID != 0) {
    logger_disable();
  } /* disable logger for MPI processes != 0 */

  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);

  /* Check if we are doing weak scaling */
  int weak = 0;
  int size = 1; // 0=small, 1=medium, 2=large, 3=very large
  if (aw != 0) {
    weak = 1;
    size = 0;
  }

  // Set up the lattice size
  if (as != 0) {
    if (strcmp(argv[as], "small") == 0) {
      size = 0;
    } else if (strcmp(argv[as], "medium") == 0) {
      size = 1;
    } else if (strcmp(argv[as], "large") == 0) {
      size = 2;
    } else if (strcmp(argv[as], "very_large") == 0) {
      size = 3;
    } else {
      error(1, 1, "read_cmdline_sombrero [sombrero.c]",
            "Arguments: [-s {small,medium,large,very_large}] [-w]");
    }
  }

  // -s: Set lattice size
  int x[4];
  if (!weak) {
    if (size == 0) {

      x[0] = 32;
      x[1] = 24;
      x[2] = 24;
      x[3] = 24;
      error(n_nodes > 1728, 1, "read_cmdline_sombrero [sombrero.c]",
            "Too many MPI ranks for a small lattice");

    } else if (size == 1) {

      x[0] = 64;
      x[1] = 48;
      x[2] = 48;
      x[3] = 48;
      error(n_nodes > 27648, 1, "read_cmdline_sombrero [sombrero.c]",
            "Too many MPI ranks for a medium lattice");

    } else if (size == 2) {

      x[0] = 96;
      x[1] = 64;
      x[2] = 64;
      x[3] = 64;
      error(n_nodes > 98304, 1, "read_cmdline_sombrero [sombrero.c]",
            "Too many MPI ranks for a large lattice");

    } else {

      x[0] = 128;
      x[1] = 96;
      x[2] = 96;
      x[3] = 96;
      error(n_nodes > 442368, 1, "read_cmdline_sombrero [sombrero.c]",
            "Too many MPI ranks for a very_large lattice");
    }
  } else {
    // Weak scaling, increase lattice size with the number of MPI ranks
    if (size == 1) {

      /* local size 24^3 = 110 592 */
      error(n_nodes < 4, 1, "main [sombrero2.c]",
            "A minimum of 4 MPI ranks is required for weak scaling with medium "
            "local volume. Try -w small or increase the number of nodes\n");
      x[0] = n_nodes;
      x[1] = 24;
      x[2] = 24;
      x[3] = 24;
      requested += 1;

    } else if (size == 2) {

      /* local size 48^3 = */
      error(n_nodes < 4, 1, "main [sombrero.c]",
            "A minimum of 4 MPI ranks is required for weak scaling with large "
            "local volume. Try -w small or increase the number of nodes\n");
      x[0] = n_nodes;
      x[1] = 48;
      x[2] = 48;
      x[3] = 48;
      requested += 1;

    } else if (size == 3) {

      /* local size 64^3 = */
      error(n_nodes < 4, 1, "main [sombrero.c]",
            "A minimum of 4 MPI ranks is required for weak scaling with "
            "very_large "
            "local volume. Try -w small or increase the number of nodes\n");
      x[0] = n_nodes;
      x[1] = 64;
      x[2] = 64;
      x[3] = 64;
      requested += 1;

    } else {

      /* local size 4^4 = 256, boundary = 240 */
      x[0] = 4;
      x[1] = 4;
      x[2] = 4;
      x[3] = 4;
      weak_scale_lattice(x, n_nodes);
    }
  }

  if (al != 0) {
    sscanf(argv[al], "%dx%dx%dx%d", x, x + 1, x + 2, x + 3);
    error(x[0] < 4 && x[1] < 4 && x[2] < 4 && x[3] < 4, 1, "main [sombrero.c]",
          "Specify lattice size as NxNxNxN, where each N is an integer larger "
          "or equal to 4.\n");
  }

  GLB_T = x[0];
  GLB_X = x[1];
  GLB_Y = x[2];
  GLB_Z = x[3];

  lprintf("GEOMETRY", 10, " Global size is %dx%dx%dx%d\n", GLB_T, GLB_X, GLB_Y,
          GLB_Z);

  /* setup process grid */
  int np[4] = {1, 1, 1, 1};

  // -p, manually partition the lattice
  if (ap != 0) {
    sscanf(argv[ap], "%dx%dx%dx%d", np, np + 1, np + 2, np + 3);
    error(np[0] < 1 && np[1] < 1 && np[2] < 1 && np[3] < 1, 1,
          "main [sombrero.c]",
          "Specify lattice size as NxNxNxN, where each N is a positive "
          "integer.\n");
  } else {
    split_lattice(np);
  }

  NP_T = np[0];
  NP_X = np[1];
  NP_Y = np[2];
  NP_Z = np[3];
  N_REP = 1;

  lprintf("GEOMETRY", 10, " Proc grid is %dx%dx%dx%d\n", NP_T, NP_X, NP_Y,
          NP_Z);
}

int main(int argc, char *argv[]) {
#if defined __SSE3__ && defined __GNUC__
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
  /* setup process communications */
  setup_process_sombrero(&argc, &argv);

  read_cmdline_sombrero(argc, argv);

  /* output version information */
  lprintf("MAIN", 0, "SOMBRERO built from HiRep commit 6c436a2\n");
  lprintf("MAIN", 0, "Base SOMBRERO commit 093dec4\n");
  lprintf("MAIN", 0,
          "SOMBRERO built using shoplifter commit 7aaec76\n");

  /* set random number variables */
  input_rlx rlx_var;
  rlx_var.rlxd_level = 1;
  rlx_var.rlxd_seed = 1;
  strcpy(rlx_var.rlxd_start, "new");
  strcpy(rlx_var.rlxd_state, "rand_state");

  /* setup lattice geometry */
  if (geometry_init() == 1) {
    finalize_process_sombrero();
    return 0;
  }

  /* setup random numbers */
  rlxd_init(rlx_var.rlxd_level,
            rlx_var.rlxd_seed +
                MPI_PID); /* use unique MPI_PID to shift seeds */

#ifdef REPR_FUNDAMENTAL
  char representation[] = "Fundamental";
#elif REPR_ADJOINT
  char representation[] = "Adjoint";
#elif REPR_SYMMETRIC
  char representation[] = "Symmetric";
#elif REPR_ANTISYMMETRIC
  char representation[] = "Antisymmetric";
#endif

  /* Predefined cases */
#if defined(GAUGE_SUN) && defined(REPR_FUNDAMENTAL) && NG == 2
  char casename[] = "Case 1";
#elif defined(GAUGE_SUN) && defined(REPR_ADJOINT) && NG == 2
  char casename[] = "Case 2";
#elif defined(GAUGE_SUN) && defined(REPR_FUNDAMENTAL) && NG == 3
  char casename[] = "Case 3";
#elif defined(GAUGE_SPN) && defined(REPR_FUNDAMENTAL) && NG == 4
  char casename[] = "Case 4";
#elif defined(GAUGE_SUN) && defined(REPR_SYMMETRIC) && NG == 3
  char casename[] = "Case 5";
#elif defined(GAUGE_SPN) && defined(REPR_ADJOINT) && NG == 4
  char casename[] = "Case 6";
#endif

  /* Initialize gauge, boundary conditions and clover  */

  /* Generate a pseudofermion field */
  float seconds;
  int real_iterations;
  {
    struct timeval start, end, etime;
    init_mc();
    spinor_field *in = alloc_spinor_field_f(1, &glat_default);
    spinor_field *out = alloc_spinor_field_f(1, &glat_default);
    flat_source(in);

    gettimeofday(&start, 0);

    int iterations = 50;
    real_iterations = cg_test(in, out, iterations);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    seconds = etime.tv_sec + etime.tv_usec / 1000000.0;
    free_spinor_field_f(out);
    /* Deallocate gauge */
    end_mc();
  }

  {
    const float Gflops = GLB_VOLUME * cg_Gflops_per_site(real_iterations);
    const float Mbytes_communicated = cg_Mbytes_communicated(real_iterations);
    const float global_GB_used =
        cg_local_GB_used(real_iterations) * NP_T * NP_X * NP_Y * NP_Z;

    lprintf("MAIN", 0, " Performed %d conjugate gradient iterations\n",
            real_iterations);
    lprintf(
        "MAIN", 0,
        " %s: %.2fe9 floating point operations and %.2fe6 bytes communicated "
        "(both directions included)\n",
        casename, Gflops, Mbytes_communicated);

    lprintf("MAIN", 0, " %s: %.2f flop per byte communicated\n", casename,
            1000 * Gflops / Mbytes_communicated);

    lprintf("MAIN", 0,
            " %s: %.2f average arithmetic intensity estimate (DRAM or L3) \n",
            casename, Gflops / global_GB_used);

    lprintf("RESULT", 0, " %s %.2f Gflops in %f seconds\n", casename, Gflops,
            seconds);
    lprintf("RESULT", 0, " %s %.2f Gflops/seconds\n", casename,
            Gflops / seconds);
  }
  /* close communications */
  finalize_process_sombrero();

  return 0;
}

/* Perform a number of iterations of the multishift conjugate gradient
 * with the Dirac operator g5Cphi_eopre_sq
 */
static int cg_test(spinor_field *in, spinor_field *out, int iterations) {

  spinor_field *k, *r, *Mk;
  spinor_field *p;
  double omega, oldomega, gamma;
  double alpha, lambda, delta;
  double innorm2;
  double *z1, *z2, *z3;
  mshift_par parameter;
  mshift_par *par = &parameter;
  double shift = 0;
  par->n = 1;
  par->shift = &shift;

  int i;
  unsigned short notconverged;

  /* fare qualche check sugli input */
#ifndef CHECK_SPINOR_MATCHING
  for (i = 0; i < par->n; ++i)
    _TWO_SPINORS_MATCHING(in, &out[i]);
#endif

  /* allocate spinors fields and aux real variables */
  p = alloc_spinor_field_f(4, in->type);
  k = p + par->n;
  r = k + 1;
  Mk = r + 1;

  z1 = malloc(sizeof(*z1));
  z2 = malloc(sizeof(*z2));
  z3 = malloc(sizeof(*z3));

  /* init recursion */
  int cgiter = 0;
  omega = 1.;
  gamma = 0.;
  innorm2 = spinor_field_sqnorm_f(in);

  /* use out[0] as initial guess */
  g5Cphi_eopre_sq(0.1, Mk, &out[0]);
  spinor_field_mul_add_assign_f(Mk, -par->shift[0], &out[0]);
  spinor_field_sub_f(r, in, Mk);

  spinor_field_copy_f(k, r);
  delta = spinor_field_sqnorm_f(r);
  z1[0] = z2[0] = 1.;
  spinor_field_copy_f(&p[0], r);

  /* cg iterations*/
  lprintf("DEBUG", 20, "%16s,%16s,%16s,%16s,%16s\n", "iteration", "alpha",
          "z3[0]_denom", "z2[0]", "delta");
  do {
    g5Cphi_eopre_sq(0.1, Mk, k);
    alpha = spinor_field_prod_re_f(k, Mk);
    if (alpha == 0)
      break;
    oldomega = omega;
    lprintf("DEBUG", 20, "%16d,", cgiter);
    lprintf("DEBUG", 20, "%16.8e,", alpha);
    omega = -delta / alpha;

    {
      double z30denom = (omega * gamma * (z1[0] - z2[0]) +
                         z1[0] * oldomega * (1. + par->shift[0] * omega));
      lprintf("DEBUG", 20, "%16.8e,", z30denom);
    }

    z3[0] = oldomega * z1[0] * z2[0] /
            (omega * gamma * (z1[0] - z2[0]) +
             z1[0] * oldomega * (1. + par->shift[0] * omega));

    lprintf("DEBUG", 20, "%16.8e,", z2[0]);
    spinor_field_mul_add_assign_f(&out[0], -omega * z3[0] / z2[0], &p[0]);

    spinor_field_mul_add_assign_f(r, omega, Mk);
    lambda = spinor_field_sqnorm_f(r);
    lprintf("DEBUG", 20, "%16.8e\n", delta);
    gamma = lambda / delta;
    delta = lambda;

    spinor_field_mul_f(k, gamma, k);
    spinor_field_add_assign_f(k, r);
    notconverged = 0; /* assume that all vectors have converged */

    /* check convergence of vectors */
    notconverged++;
    spinor_field_mul_f(&p[0], gamma * z3[0] * z3[0] / (z2[0] * z2[0]), &p[0]);
    spinor_field_mul_add_assign_f(&p[0], z3[0], r);
    z1[0] = z2[0];
    z2[0] = z3[0];

    ++cgiter;
  } while (cgiter < iterations);

  {
    double output_norm = spinor_field_sqnorm_f(Mk);
    g5Cphi_eopre_sq(0.1, Mk, &out[0]);
    spinor_field_mul_add_assign_f(Mk, -par->shift[0], &out[0]);
    spinor_field_sub_f(Mk, Mk, in);

    double input_norm = spinor_field_sqnorm_f(in);
    double residual_norm = spinor_field_sqnorm_f(Mk);
    double ratio = residual_norm / input_norm;
    lprintf("RESULT", 20, "Deviation from expected %16.8e\n", ratio);
    lprintf(
        "RESULT", 20,
        "Residual norm: %16.8e, norm of input: %16.8e, output norm: %16.8e\n",
        residual_norm, input_norm, output_norm); // DEBUG
    error(ratio > 1e-8, 1, "main [sombrero.c]",
          "Result of CG inversion incorrect\n");
  }

  /* free memory */
  free_spinor_field_f(p);
  free(z1);
  free(z2);
  free(z3);

  /* return number of cg iter */
  return cgiter;
}
