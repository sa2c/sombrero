/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "dirac.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include "communications.h"
#include <math.h>

void gaussian_spinor_field(spinor_field *s) {
	const double c1=1./sqrt(2.);
	int i;
	geometry_descriptor *type = s->type;
	for(i=0;i<type->local_master_pieces;i++)
 	  gauss((double*)(s->ptr+(type->master_start[i]-type->master_shift)),(type->master_end[i]-type->master_start[i]+1)*sizeof(suNf_spinor)/sizeof(double));
	spinor_field_mul_f(s,c1,s);
	apply_BCs_on_spinor_field(s);
}

void flat_source(spinor_field *s) {
	_ONE_SPINOR_FOR(s) {
    	_spinor_one_f(*_SPINOR_PTR(s));
    }
	//if(CID == 0){
	// 	_FIELD_AT(s,0)->c[0].c[0].re = 1.0;
	//}
	apply_BCs_on_spinor_field(s);

	start_sf_sendrecv(s);
	complete_sf_sendrecv(s);
}

void gaussian_spinor_field_flt(spinor_field_flt *s) {
	const float c1=(float)(1./sqrt(2.));
	int i;
	geometry_descriptor *type = s->type;
	for(i=0;i<type->local_master_pieces;i++)
 	  gauss_flt((float*)(s->ptr+(type->master_start[i]-type->master_shift)),(type->master_end[i]-type->master_start[i]+1)*sizeof(suNf_spinor_flt)/sizeof(float));
	spinor_field_mul_f_flt(s,c1,s);
	apply_BCs_on_spinor_field_flt(s);
}


