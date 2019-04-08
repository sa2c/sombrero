/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "complex.h"
#include "spinor_field.h"


typedef struct _mshift_par {
   int n; /* number of shifts */
   double *shift;
   double err2; /* relative error of the solutions */
   int max_iter; /* maximum number of iterations: 0 => infinity */
	 void *add_par; /* additional parameters for specific inverters */
} mshift_par;


#endif
