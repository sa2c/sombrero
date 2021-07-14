/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random.h
* 
* Pseudorandom number, matrices and fields
*
*******************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include "suN.h"
#include "spinor_field.h"

#include "ranlux.h" /* included here for convenience */


void random_suNg_unit_vector(suNg_vector *v);

void random_suNg(suNg *u);

void random_u(suNg_field *gf);



#endif
