/*******************************************************************************
*
* File suN_types.h
*
* Type definitions for SU(N) matrices and spinors
*
*******************************************************************************/

#ifndef SUN_TYPES_H
#define SUN_TYPES_H

#include "complex.h"

/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

typedef struct
{
   complex c[3];
} suNg_vector;

typedef struct
{
   complex_flt c[3];
} suNg_vector_flt;

typedef struct
{
   complex c[9];
} suNg;

typedef struct
{
   complex_flt c[9];
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
   double c[8];
} suNg_algebra_vector;

typedef struct
{
   float c[8];
} suNg_algebra_vector_flt;

#define NF 6
/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

typedef struct
{
   complex c[6];
} suNf_vector;

typedef struct
{
   complex_flt c[6];
} suNf_vector_flt;

typedef struct
{
   complex c[36];
} suNf;

typedef struct
{
   complex_flt c[36];
} suNf_flt;

typedef suNf suNfc;

typedef suNf_flt suNfc_flt;

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
