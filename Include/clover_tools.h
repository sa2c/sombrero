/***************************************************************************\
* Copyright (c) 2016, Martin Hansen                                         *
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef CLOVER_TOOLS_H
#define CLOVER_TOOLS_H

#include "suN.h"
#include "spinor_field.h"

#if (defined(WITH_CLOVER) && defined(UPDATE_EO))
#define WITH_CLOVER_EO
#endif

void compute_ldl_decomp(double);
void compute_clover_term();
void clover_init(double);

#endif
