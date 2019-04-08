/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File Dphi.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given
* double-precision spinor field
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "error.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "spinor_field.h"
#include "geometry.h"
#include "communications.h"
#include "memory.h"
#include "clover_tools.h"

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */


/*
 * Init of Dphi
 */

int init_dirac=1;
spinor_field *gtmp=NULL;
spinor_field *etmp=NULL;
spinor_field *otmp=NULL;
spinor_field *otmp2=NULL;

static void free_mem() {
    if (gtmp!=NULL) { free_spinor_field_f(gtmp); etmp=NULL; }
    if (etmp!=NULL) { free_spinor_field_f(etmp); etmp=NULL; }
    if (otmp!=NULL) { free_spinor_field_f(otmp); otmp=NULL; }
    if (otmp2!=NULL) { free_spinor_field_f(otmp2); otmp2=NULL; }
    init_dirac=1;
}

void init_Dirac() {
    if (init_dirac) {
        gtmp=alloc_spinor_field_f(1,&glattice);
        etmp=alloc_spinor_field_f(1,&glat_even);
        otmp=alloc_spinor_field_f(1,&glat_odd);
        otmp2=alloc_spinor_field_f(1,&glat_odd);
        atexit(&free_mem);
        init_dirac=0;
    }
}


/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	MVMcounter=0; /* reset counter */

	return res;
}


/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r,u,s)\
    _suNf_multiply(vtmp,(u),(s));\
    _vector_mulc_f((r),eitheta[0],vtmp)

#define _suNf_theta_T_inverse_multiply(r,u,s)\
    _suNf_inverse_multiply(vtmp,(u),(s));\
    _vector_mulc_star_f((r),eitheta[0],vtmp)

#else

#define _suNf_theta_T_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_T_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[1],vtmp)

#define _suNf_theta_X_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[1],vtmp)

#else

#define _suNf_theta_X_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_X_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[2],vtmp)

#define _suNf_theta_Y_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[2],vtmp)

#else

#define _suNf_theta_Y_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Y_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[3],vtmp)

#define _suNf_theta_Z_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[3],vtmp)

#else

#define _suNf_theta_Z_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Z_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif




/*
 * This function defines the massless Dirac operator
 * It can act on spinors defined on the whole lattice
 * or on spinors with definite parity
 *
 * Flop count = VOL*(16*MMUL + 76*NF), where
 */
void Dphi_(spinor_field *out, spinor_field *in)
{

   error((in==NULL)||(out==NULL),1,"Dphi_ [Dphi.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"Dphi_ [Dphi.c]",
         "Input and output fields must be different");

#ifndef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd,1,"Dphi_ [Dphi.c]", "Spinors don't match! (1)");
   error(out->type==&glat_odd && in->type!=&glat_even,1,"Dphi_ [Dphi.c]", "Spinors don't match! (2)");
   error(out->type==&glattice && in->type!=&glattice,1,"Dphi_ [Dphi.c]", "Spinors don't match! (3)");
#endif


   ++MVMcounter; /* count matrix calls */
   if(out->type==&glattice) ++MVMcounter;

/************************ loop over all lattice sites *************************/
   /* start communication of input spinor field */
   start_sf_sendrecv(in);

  _PIECE_FOR(out->type,ixp) {
     if(ixp==out->type->inner_master_pieces) {
       _OMP_PRAGMA( master )
       /* wait for spinor to be transfered */
       complete_sf_sendrecv(in);
       _OMP_PRAGMA( barrier )
     }
     _SITE_FOR(out->type,ixp,ix) {

       int iy;
       suNf *up,*um;
       suNf_vector psi,chi;
       suNf_spinor *r=0,*sp,*sm;
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
       suNf_vector vtmp;
#endif

       r=_FIELD_AT(out,ix);

       /******************************* direction +0 *********************************/

       iy=iup(ix,0);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,0);

       _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[0]=chi;
       (*r).c[2]=chi;

       _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[1]=chi;
       (*r).c[3]=chi;

       /******************************* direction -0 *********************************/

       iy=idn(ix,0);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,0);

       _vector_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       _vector_sub_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       /******************************* direction +1 *********************************/

       iy=iup(ix,1);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,1);

       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);

       _vector_i_add_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       /******************************* direction -1 *********************************/

       iy=idn(ix,1);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,1);

       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       _vector_i_sub_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       /******************************* direction +2 *********************************/

       iy=iup(ix,2);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,2);

       _vector_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_add_assign_f((*r).c[3],chi);

       _vector_sub_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       /******************************* direction -2 *********************************/

       iy=idn(ix,2);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,2);

       _vector_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       _vector_add_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_add_assign_f((*r).c[2],chi);

       /******************************* direction +3 *********************************/

       iy=iup(ix,3);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,3);

       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       _vector_i_sub_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       /******************************* direction -3 *********************************/

       iy=idn(ix,3);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,3);

       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);

       /******************************** end of loop *********************************/

       _spinor_mul_f(*r,-0.5,*r);

     } /* SITE_FOR */
   } /* PIECE FOR */
}


#ifdef WITH_CLOVER

/*************************************************
 * Dirac operators with clover term:             *
 * Cphi = Dphi + clover                          *
 * Cphi_eopre = D_ee - D_eo D_oo^-1 D_oe         *
 * Cphi_diag = D_oo or D_ee                      *
 * Cphi_diag_inv = D_oo^-1 or D_ee^-1            *
 *************************************************/

#if defined(REPR_ADJOINT)
typedef suNfc clover_type;
#define clover_inverse_multiply(a,b,c) _suNfc_inverse_multiply(a,b,c) 
#define clover_multiply(a,b,c) _suNfc_multiply(a,b,c) 
#else
typedef suNffull clover_type;
#define clover_inverse_multiply(a,b,c) _suNffull_inverse_multiply(a,b,c) 
#define clover_multiply(a,b,c) _suNffull_multiply(a,b,c) 
#endif

/* Flop count = VOL*(16*NF*NF + 24*NF), where
 */
static void Cphi_(double mass, spinor_field *dptr, spinor_field *sptr, int assign)
{
	// Correct mass term
	mass = (4.+mass);

	// Loop over local sites
	_MASTER_FOR(dptr->type,ix)
	{
		suNf_vector v1, v2;
		suNf_spinor *out, *in, tmp;
		clover_type *s0, *s1, *s2, *s3;

		// Field pointers
		out = _FIELD_AT(dptr,ix);
		in = _FIELD_AT(sptr,ix);
		s0 = _4FIELD_AT(cl_term,ix,0);
		s1 = _4FIELD_AT(cl_term,ix,1);
		s2 = _4FIELD_AT(cl_term,ix,2);
		s3 = _4FIELD_AT(cl_term,ix,3);

		// Component 0
		clover_multiply(v1, *s0, in->c[0]);
		clover_multiply(v2, *s1, in->c[1]);
		_vector_add_f(tmp.c[0], v1, v2);

		// Component 1
		clover_inverse_multiply(v1, *s1, in->c[0]);
		clover_multiply(v2, *s0, in->c[1]);
		_vector_sub_f(tmp.c[1], v1, v2);

		// Component 2
		clover_multiply(v1, *s2, in->c[2]);
		clover_multiply(v2, *s3, in->c[3]);
		_vector_add_f(tmp.c[2], v1, v2);

		// Component 3
		clover_inverse_multiply(v1, *s3, in->c[2]);
		clover_multiply(v2, *s2, in->c[3]);
		_vector_sub_f(tmp.c[3], v1, v2);

		// Add mass
		_spinor_mul_add_assign_f(tmp, mass, *in);

		// Store
		if(assign)
		{
			_spinor_add_assign_f(*out, tmp);
		}
		else
		{
			*out = tmp;
		}
	}
}

/* Flop count = VOL*(56*NF*NF + 32*NF), where
 */
static void Cphi_inv_(double mass, spinor_field *dptr, spinor_field *sptr, int assign)
{
	int N = 2*NF;
	mass = (4.+mass);

	// Update LDL decomposition
	compute_ldl_decomp(mass);

	// Loop over local sites
	_MASTER_FOR(dptr->type,ix)
	{
		complex *up, *dn, *x, c;
		suNf_spinor *out, *in, tmp;
		int n;

		// Field pointers
		up = _FIELD_AT(cl_ldl,ix)->up;
		dn = _FIELD_AT(cl_ldl,ix)->dn;
		out = _FIELD_AT(dptr,ix);
		in = _FIELD_AT(sptr,ix);

		// tmp = in
		tmp = *in;
		x = (complex*)&tmp;

		// Forward substitution
		for(int i = 0; i < N; i++)
		{
			for(int k = 0; k < i; k++)
			{
				n = i*(i+1)/2+k;
				_complex_mul_sub_assign(x[i], up[n], x[k]);
				_complex_mul_sub_assign(x[i+N], dn[n], x[k+N]);
			}
		}

		// Backward substitution
		for(int i = N-1; i >= 0; i--)
		{
			n = i*(i+1)/2+i;
			_complex_mulr(x[i], 1./up[n].re, x[i]);
			_complex_mulr(x[i+N], 1./dn[n].re, x[i+N]);
			for(int k = i+1; k < N; k++)
			{
				n = k*(k+1)/2+i;

				c.re = up[n].re;
				c.im = -up[n].im;
				_complex_mul_sub_assign(x[i], c, x[k]);

				c.re = dn[n].re;
				c.im = -dn[n].im;
				_complex_mul_sub_assign(x[i+N], c, x[k+N]);
			}
		}

		// Store
		if(assign)
		{
			_spinor_add_assign_f(*out, tmp);
		}
		else
		{
			*out = tmp;
		}
	}
}

/* Flop count = VOL*(16*MMUL + 72NF*NF + 128*NF), where
 * MMUL=4*NF*NF - 2*NF for adjoint, 8*NF*NF - 2*NF otherwise
 */
void Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr)
{
	if(init_dirac)
	{
		init_Dirac();
	}

	apply_BCs_on_spinor_field(sptr);
	Dphi_(otmp, sptr);
	Cphi_inv_(mass, otmp, otmp, 0);
	apply_BCs_on_spinor_field(otmp);
	Dphi_(dptr, otmp);
	spinor_field_minus_f(dptr, dptr);
	Cphi_(mass, dptr, sptr, 1);
	apply_BCs_on_spinor_field(dptr);
}

void g5Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr)
{
	Cphi_eopre(mass, dptr, sptr);
	spinor_field_g5_assign_f(dptr);
}

/* Flop count = VOL*(32*MMUL + 144*NF*NF + 256*NF), where
 * MMUL=4*NF*NF - 2*NF for adjoint, 8*NF*NF - 2*NF otherwise
 */
void g5Cphi_eopre_sq(double mass, spinor_field *dptr, spinor_field *sptr)
{
	if(init_dirac)
	{
		init_Dirac();
	}

	g5Cphi_eopre(mass, etmp, sptr);
	g5Cphi_eopre(mass, dptr, etmp);
}

#endif //#ifdef WITH_CLOVER
