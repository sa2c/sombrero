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

#if !defined(WITH_QUATERNIONS) && !defined(GAUGE_SON)
/* generates a random SU(N) matrix via SU(2) rotations */
static void rotate(suNg_vector *pu1, suNg_vector *pu2, double s[4]) /* same as in cabmar */
{
    double complex z1,z2;
    double complex *cu1, *cu2;

    cu1 = &((*pu1).c[0]);
    cu2 = &((*pu2).c[0]);

    for (int i=0; i<NG; ++i) {
        z1=s[0]*creal(*cu1)-s[1]*cimag(*cu2)+s[2]*creal(*cu2)-s[3]*cimag(*cu1)+I*(s[0]*cimag(*cu1)+s[1]*creal(*cu2)+s[2]*cimag(*cu2)+s[3]*creal(*cu1));
        z2=s[0]*creal(*cu2)-s[1]*cimag(*cu1)-s[2]*creal(*cu1)+s[3]*cimag(*cu2)+I*(s[0]*cimag(*cu2)+s[1]*creal(*cu1)-s[2]*cimag(*cu1)-s[3]*creal(*cu2));
        (*cu1) = z1;
        (*cu2) = z2;
        ++cu1;
        ++cu2;
    }
}
#endif 

#ifndef GAUGE_SON
#ifdef GAUGE_SPN
#define _VARNAME r
#else 
#define _VARNAME u
#endif
void random_suNg(suNg *_VARNAME) {
#undef _VARNAME
#ifdef GAUGE_SPN
    suNgfull ut, *u;
    _suNg_expand(ut,*r);
    u=&ut;
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
        for (int i=0; i<NG*NG/2; ++i) { r->c[i]=ut.c[i]; }
#endif
}
#else // vv GAUGE_SON case vv
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

void random_suNf(suNf *u) {
    suNf tmp;
    double gr[NF*NF];
    int i;
    do {
        gauss(gr,NF*NF);
        for (i=0;i<NF*NF;i++){
            tmp.c[i]=gr[i];
        }
    } while (!project_to_suNg_real(u,&tmp));
}


#endif

