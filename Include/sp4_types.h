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
   complex c[4];
} suNg_vector;

typedef struct
{
   complex_flt c[4];
} suNg_vector_flt;

typedef struct
{
   complex c[8];
} suNg;

typedef struct
{
   complex_flt c[8];
} suNg_flt;

typedef suNg suNgc;

typedef suNg_flt suNgc_flt;

typedef struct
{
   complex c[16];
} suNgfull;

typedef struct
{
   complex_flt c[16];
} suNgfull_flt;

typedef suNgfull suNg_FMAT;

typedef suNgfull_flt suNg_FMAT_flt;

typedef suNgfull suNgfullc;

typedef suNgfull_flt suNgfullc_flt;

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
   double c[10];
} suNg_algebra_vector;

typedef struct
{
   float c[10];
} suNg_algebra_vector_flt;

#define NF 4
/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

typedef struct
{
   complex c[4];
} suNf_vector;

typedef struct
{
   complex_flt c[4];
} suNf_vector_flt;

typedef struct
{
   complex c[8];
} suNf;

typedef struct
{
   complex_flt c[8];
} suNf_flt;

typedef suNf suNfc;

typedef suNf_flt suNfc_flt;

typedef struct
{
   complex c[16];
} suNffull;

typedef struct
{
   complex_flt c[16];
} suNffull_flt;

typedef suNffull suNf_FMAT;

typedef suNffull_flt suNf_FMAT_flt;

typedef suNffull suNffullc;

typedef suNffull_flt suNffullc_flt;

typedef suNffull suNffullfull;

typedef struct
{
   suNf_vector c[4];
} suNf_spinor;

typedef struct
{
   suNf_vector_flt c[4];
} suNf_spinor_flt;


#endif
