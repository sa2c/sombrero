/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef IO_H
#define IO_H
#include "input_par.h"
#include "spinor_field.h"
#include <stdio.h>
#include "suN.h"
#include "update.h"




void read_gauge_field_nocheck(char filename[]);

#ifdef GAUGE_SPN
#endif




void read_spinor_field_ascii(char filename[],spinor_field * sf);

void print_mat(suNg* mat);
#endif
