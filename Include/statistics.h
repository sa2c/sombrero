/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File statistics.h
* 
* Functions for statistical analysis of data series
*
*******************************************************************************/

#ifndef EXTRAS_H
#define EXTRAS_H

double average(int n,double a[]);
double sigma0(int n,double a[]);
void auto_corr(int n,double a[],int tmax,double gamma[]);
double sigma(int n,double a[],double *tau,int *flag);



#endif

 
