/* This is part of the ranlux 3.3 distribution by Martin Luesher */

/*******************************************************************************
*
* File ranlxd.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Random number generator "ranlxd". See the notes 
*
*   "User's guide for ranlxs and ranlxd v3.2" (December 2005)
*
*   "Algorithms used in ranlux v3.0" (May 2001)
*
* for a detailed description
*
* The externally accessible functions are 
*
*   void ranlxd(double r[],int n)
*     Computes the next n double-precision random numbers and 
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
* 
*   void rlxd_init(int level,int seed)
*     Initialization of the generator
*
*   int rlxd_size(void)
*     Returns the number of integers required to save the state of
*     the generator
*
*   void rlxd_get(int state[])
*     Extracts the current state of the generator and stores the 
*     information in the array state[N] where N>=rlxd_size()
*
*   void rlxd_reset(int state[])
*     Resets the generator to the state defined by the array state[N]
*
*******************************************************************************/

#define RANLXD_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlux.h"

#if (defined SSE)

typedef struct
{
      float c1, c2, c3, c4;
} vec_t __attribute__((aligned(16)));

typedef struct
{
      vec_t c1, c2;
} dble_vec_t __attribute__((aligned(16)));

static int init = 0, pr, prm, ir, jr, is, is_old, next[96];
static vec_t one, one_bit, carry;

static union {
      dble_vec_t vec[12];
      float num[96];
} x __attribute__((aligned(16)));

#define STEP(pi, pj)                                         \
      __asm__ __volatile__("movaps %4, %%xmm4 \n\t"          \
                           "movaps %%xmm2, %%xmm3 \n\t"      \
                           "subps %2, %%xmm4 \n\t"           \
                           "movaps %%xmm1, %%xmm5 \n\t"      \
                           "cmpps $0x6, %%xmm4, %%xmm2 \n\t" \
                           "andps %%xmm2, %%xmm5 \n\t"       \
                           "subps %%xmm3, %%xmm4 \n\t"       \
                           "andps %%xmm0, %%xmm2 \n\t"       \
                           "addps %%xmm4, %%xmm5 \n\t"       \
                           "movaps %%xmm5, %0 \n\t"          \
                           "movaps %5, %%xmm6 \n\t"          \
                           "movaps %%xmm2, %%xmm3 \n\t"      \
                           "subps %3, %%xmm6 \n\t"           \
                           "movaps %%xmm1, %%xmm7 \n\t"      \
                           "cmpps $0x6, %%xmm6, %%xmm2 \n\t" \
                           "andps %%xmm2, %%xmm7 \n\t"       \
                           "subps %%xmm3, %%xmm6 \n\t"       \
                           "andps %%xmm0, %%xmm2 \n\t"       \
                           "addps %%xmm6, %%xmm7 \n\t"       \
                           "movaps %%xmm7, %1"               \
                           : "=m"((*pi).c1),                 \
                             "=m"((*pi).c2)                  \
                           : "m"((*pi).c1),                  \
                             "m"((*pi).c2),                  \
                             "m"((*pj).c1),                  \
                             "m"((*pj).c2)                   \
                           : "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7")

static void error(int no)
{
      switch (no)
      {
      case 1:
            printf("Error in subroutine rlxd_init\n");
            printf("Bad choice of luxury level (should be 1 or 2)\n");
            break;
      case 2:
            printf("Error in subroutine rlxd_init\n");
            printf("Bad choice of seed (should be between 1 and 2^31-1)\n");
            break;
      case 3:
            printf("Error in rlxd_get\n");
            printf("Undefined state (ranlxd is not initialized\n");
            break;
      case 5:
            printf("Error in rlxd_reset\n");
            printf("Unexpected input data\n");
            break;
      }
      printf("Program aborted\n");
      exit(0);
}


static void define_constants(void)
{
      int k;
      float b;

      one.c1 = 1.0f;
      one.c2 = 1.0f;
      one.c3 = 1.0f;
      one.c4 = 1.0f;

      b = (float)(ldexp(1.0, -24));
      one_bit.c1 = b;
      one_bit.c2 = b;
      one_bit.c3 = b;
      one_bit.c4 = b;

      for (k = 0; k < 96; k++)
      {
            next[k] = (k + 1) % 96;
            if ((k % 4) == 3)
                  next[k] = (k + 5) % 96;
      }
}

void rlxd_init(int level, int seed)
{
      int i, k, l;
      int ibit, jbit, xbit[31];
      int ix, iy;

      define_constants();

      if (level == 1)
            pr = 202;
      else if (level == 2)
            pr = 397;
      else
            error(1);

      i = seed;

      for (k = 0; k < 31; k++)
      {
            xbit[k] = i % 2;
            i /= 2;
      }

      if ((seed <= 0) || (i != 0))
            error(2);

      ibit = 0;
      jbit = 18;

      for (i = 0; i < 4; i++)
      {
            for (k = 0; k < 24; k++)
            {
                  ix = 0;

                  for (l = 0; l < 24; l++)
                  {
                        iy = xbit[ibit];
                        ix = 2 * ix + iy;

                        xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                        ibit = (ibit + 1) % 31;
                        jbit = (jbit + 1) % 31;
                  }

                  if ((k % 4) != i)
                        ix = 16777215 - ix;

                  x.num[4 * k + i] = (float)(ldexp((double)(ix), -24));
            }
      }

      carry.c1 = 0.0f;
      carry.c2 = 0.0f;
      carry.c3 = 0.0f;
      carry.c4 = 0.0f;

      ir = 0;
      jr = 7;
      is = 91;
      is_old = 0;
      prm = pr % 12;
      init = 1;
}





#else

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
      int c1, c2, c3, c4;
} vec_t;

typedef struct
{
      vec_t c1, c2;
} dble_vec_t;

static int init = 0, pr, prm, ir, jr, is, is_old, next[96];
static double one_bit;
static vec_t carry;

static union {
      dble_vec_t vec[12];
      int num[96];
} x;

#define STEP(pi, pj)                            \
      d = (*pj).c1.c1 - (*pi).c1.c1 - carry.c1; \
      (*pi).c2.c1 += (d < 0);                   \
      d += BASE;                                \
      (*pi).c1.c1 = d & MASK;                   \
      d = (*pj).c1.c2 - (*pi).c1.c2 - carry.c2; \
      (*pi).c2.c2 += (d < 0);                   \
      d += BASE;                                \
      (*pi).c1.c2 = d & MASK;                   \
      d = (*pj).c1.c3 - (*pi).c1.c3 - carry.c3; \
      (*pi).c2.c3 += (d < 0);                   \
      d += BASE;                                \
      (*pi).c1.c3 = d & MASK;                   \
      d = (*pj).c1.c4 - (*pi).c1.c4 - carry.c4; \
      (*pi).c2.c4 += (d < 0);                   \
      d += BASE;                                \
      (*pi).c1.c4 = d & MASK;                   \
      d = (*pj).c2.c1 - (*pi).c2.c1;            \
      carry.c1 = (d < 0);                       \
      d += BASE;                                \
      (*pi).c2.c1 = d & MASK;                   \
      d = (*pj).c2.c2 - (*pi).c2.c2;            \
      carry.c2 = (d < 0);                       \
      d += BASE;                                \
      (*pi).c2.c2 = d & MASK;                   \
      d = (*pj).c2.c3 - (*pi).c2.c3;            \
      carry.c3 = (d < 0);                       \
      d += BASE;                                \
      (*pi).c2.c3 = d & MASK;                   \
      d = (*pj).c2.c4 - (*pi).c2.c4;            \
      carry.c4 = (d < 0);                       \
      d += BASE;                                \
      (*pi).c2.c4 = d & MASK

static void error(int no)
{
      switch (no)
      {
      case 0:
            printf("Error in rlxd_init\n");
            printf("Arithmetic on this machine is not suitable for ranlxd\n");
            break;
      case 1:
            printf("Error in subroutine rlxd_init\n");
            printf("Bad choice of luxury level (should be 1 or 2)\n");
            break;
      case 2:
            printf("Error in subroutine rlxd_init\n");
            printf("Bad choice of seed (should be between 1 and 2^31-1)\n");
            break;
      case 3:
            printf("Error in rlxd_get\n");
            printf("Undefined state (ranlxd is not initialized)\n");
            break;
      case 4:
            printf("Error in rlxd_reset\n");
            printf("Arithmetic on this machine is not suitable for ranlxd\n");
            break;
      case 5:
            printf("Error in rlxd_reset\n");
            printf("Unexpected input data\n");
            break;
      }
      printf("Program aborted\n");
      exit(0);
}


static void define_constants(void)
{
      int k;

      one_bit = ldexp(1.0, -24);

      for (k = 0; k < 96; k++)
      {
            next[k] = (k + 1) % 96;
            if ((k % 4) == 3)
                  next[k] = (k + 5) % 96;
      }
}

void rlxd_init(int level, int seed)
{
      int i, k, l;
      int ibit, jbit, xbit[31];
      int ix, iy;

      if ((INT_MAX < 2147483647) || (FLT_RADIX != 2) || (FLT_MANT_DIG < 24) ||
          (DBL_MANT_DIG < 48))
            error(0);

      define_constants();

      if (level == 1)
            pr = 202;
      else if (level == 2)
            pr = 397;
      else
            error(1);

      i = seed;

      for (k = 0; k < 31; k++)
      {
            xbit[k] = i % 2;
            i /= 2;
      }

      if ((seed <= 0) || (i != 0))
            error(2);

      ibit = 0;
      jbit = 18;

      for (i = 0; i < 4; i++)
      {
            for (k = 0; k < 24; k++)
            {
                  ix = 0;

                  for (l = 0; l < 24; l++)
                  {
                        iy = xbit[ibit];
                        ix = 2 * ix + iy;

                        xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                        ibit = (ibit + 1) % 31;
                        jbit = (jbit + 1) % 31;
                  }

                  if ((k % 4) != i)
                        ix = 16777215 - ix;

                  x.num[4 * k + i] = ix;
            }
      }

      carry.c1 = 0;
      carry.c2 = 0;
      carry.c3 = 0;
      carry.c4 = 0;

      ir = 0;
      jr = 7;
      is = 91;
      is_old = 0;
      prm = pr % 12;
      init = 1;
}





#endif
