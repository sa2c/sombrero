/*******************************************************************************
*
* File suN_types.h
*
* Type definitions for SU(N) matrices and spinors
*
*******************************************************************************/

#ifndef SUN_TYPES_H
#define SUN_TYPES_H

#include "hr_complex.h"

#define NG 2
/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

typedef struct
{
   double complex c[2];
} suNg_vector;

typedef struct
{
   float complex c[2];
} suNg_vector_flt;

typedef struct
{
   double complex c[4];
} suNg;

typedef struct
{
   float complex c[4];
} suNg_flt;

typedef suNg suNgc;

typedef suNg_flt suNgc_flt;

typedef suNg suNg_FMAT;

typedef suNg_flt suNg_FMAT_flt;

typedef struct
{
   suNg_vector c[4];
} suNg_spinor;

typedef struct
{
   suNg_vector_flt c[4];
} suNg_spinor_flt;

typedef struct
{
   double c[3];
} suNg_algebra_vector;

typedef struct
{
   float c[3];
} suNg_algebra_vector_flt;

#define NF 3
/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

typedef struct
{
   double complex c[3];
} suNf_vector;

typedef struct
{
   float complex c[3];
} suNf_vector_flt;

typedef struct
{
   double complex c[9];
} suNfc;

typedef struct
{
   float complex c[9];
} suNfc_flt;

typedef struct
{
   double c[9];
} suNf;

typedef struct
{
   float c[9];
} suNf_flt;

typedef suNf suNf_FMAT;

typedef suNf_flt suNf_FMAT_flt;

typedef suNf suNffull;

typedef struct
{
   suNf_vector c[4];
} suNf_spinor;

typedef struct
{
   suNf_vector_flt c[4];
} suNf_spinor_flt;


#endif
