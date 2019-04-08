/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random_suNg.c
*
* Generation of uniformly distributed SU(Ng) matrices
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "utils.h"

extern void random_su2(double rho,double s[]);

void gaussian_suNg_vector(suNg_vector *v)
{
   gauss((double*)v, sizeof(suNg_vector)/sizeof(double));
}

/* generates a random SU(N) matrix via SU(2) rotations */
static void rotate(suNg_vector *pu1, suNg_vector *pu2, double s[4]) /* same as in cabmar */
{
	  complex z1,z2;
	  complex *cu1, *cu2;
  
	  cu1 = &((*pu1).c[0]);
	  cu2 = &((*pu2).c[0]);
						  
	  for (int i=0; i<NG; ++i) {
      z1.re=s[0]*(*cu1).re-s[1]*(*cu2).im+s[2]*(*cu2).re-s[3]*(*cu1).im;
      z1.im=s[0]*(*cu1).im+s[1]*(*cu2).re+s[2]*(*cu2).im+s[3]*(*cu1).re;
      z2.re=s[0]*(*cu2).re-s[1]*(*cu1).im-s[2]*(*cu1).re+s[3]*(*cu2).im;
      z2.im=s[0]*(*cu2).im+s[1]*(*cu1).re-s[2]*(*cu1).im-s[3]*(*cu2).re;
      (*cu1) = z1;
      (*cu2) = z2;
      ++cu1;
      ++cu2;
	  }
}

#ifndef GAUGE_SON
#ifdef GAUGE_SPN
void random_suNg(suNg *r) {
  suNgfull ut, *u;
  _suNg_expand(ut,*r);
  u=&ut;
#else
void random_suNg(suNg *u) {
#endif
#ifdef WITH_QUATERNIONS
  random_su2(0.,u->c);
#else
  double s[4];
  suNg_vector *pu1=(suNg_vector*)(u);
	
#ifdef GAUGE_SPN
  _suNgfull_unit(*u);
#else
  _suNg_unit(*u);
#endif
  
  for (int i=0; i<NG; ++i) {
    suNg_vector *pu2 = pu1 + 1;
    for (int j=i+1; j<NG; ++j) {
  #ifdef GAUGE_SPN
  	if( (i < NG/2 && j < NG/2) ){
		random_su2(0.0,s);
  		rotate(pu1, pu2, s);
  		s[3] *=-1.;
  		s[1] *=-1.;
  		rotate(pu1+NG/2, pu2+NG/2,s);
  	} else if( ( i >= NG/2 && j > NG/2) ){
        random_su2(0.0,s);
  		rotate(pu1,pu2,s);
  		s[3] *=-1.;		
  		s[1] *=-1.;
  		rotate(pu1-NG/2,pu2-NG/2,s);
  	}
    else if( ( j == i + NG/2) ){
        random_su2(0.0,s);
        rotate(pu1,pu2,s) ;
    }
  	else ;
  #else
    random_su2(0.0,s);
  	rotate(pu1,pu2,s);
  #endif
      ++pu2; 
    } 
	  ++pu1; 
  }
#endif //WITH_QUATERNIONS
#ifdef GAUGE_SPN
	for (int i=0; i<NG*NG/2; ++i) { r->c[i].re=ut.c[i].re; r->c[i].im=ut.c[i].im; }
#endif
}

#else // GAUGE_SON case
void random_suNg(suNg *u) {
  suNg tmp;
  double gr[NG*NG];
  int i;
  do {
    gauss(gr,NG*NG);
    for (i=0;i<NG*NG;i++){
      tmp.c[i]=gr[i];
    }
  } while (!project_to_suNg_real(u,&tmp));
}

#endif
