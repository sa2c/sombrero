/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"

void project_gauge_field(void);
void random_su2(double rho,double s[]);

void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);

void flat_source(spinor_field *s);

#endif
