/*******************************************************************
   Copyright (C) 2001-2 Leo Breiman, Adele Cutler, Andy Liaw and Mathew Wiener
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.                            
*******************************************************************/

/* Template of Fortran subroutines to be called from the C wrapper */
extern void F77_NAME(createclass)(double *, int *, int *, int *, int *,
				  double *, double *, double *, int *, int *);
extern void F77_NAME(zerm)(int *,     int *, int *);
extern void F77_NAME(zermd)(double *, int *, int *);
extern void F77_NAME(zerv)(int *,     int *);
extern void F77_NAME(zervr)(double *, int *);
extern void F77_NAME(prep)(int *, int *, int *, int *, double *, double *, 
			   int*, double *); 
extern void F77_NAME(makea)(double *, int *, int *, int *, int *, double *,
			    int *, int *, int *);
extern void F77_NAME(eqm)(int *, int *, int *, int *);
extern void F77_NAME(moda)(int *, int *, int *, int *, int *, int *, int *,
			   int *, int *);
extern void F77_NAME(buildtree)(int *, int *, int *, int *, int *, int *, int
				*, int *, int *, int *, int *, double *, int
				*, int *, int *, double *, double *, double
				*, int *, int *, int *, int *, int *, int *,
				int *, int *, int *, int *, int *, double *,
				double *, double *, double *, int *, int *,
				int *); 
extern void F77_NAME(xtranslate)(double *, int *, int *, int *, int *, int *,
				 int *, double *, int *, int *, int *); 
extern void F77_NAME(testreebag)(double *, int *, int *, int *, int *,
				 double *, int *, int *, int *, int *, int *,
				 int *, int *, int *, int *, int *); 
extern void F77_NAME(oob)(int *, int *, int *, int *, int *, int *, int *,
			  int *, double *, double *, double *, double *, int
			  *, double *);
extern void F77_NAME(permobmr)(int *, double *, double *, double *, int *,
			       int *, int *);
extern void F77_NAME(locateout)(double *, int *, int *, int *, int *c, int *,
				int *,	double *, double *, int *, int *);
extern void F77_NAME(finishimp)(double *, int *, int *, int * , int *, int *,
				int *, double *, double *, double *, double *, 
				double *, int *, int *, double *);

/*************************************************************
  Subroutines for the regression RF
************************************************************/

extern void F77_NAME(rbuildtree)(double *xb, double *yb, double *yl, int *mdim,
				int *nls, int *nsample, int *treemap, 
				int *jdex, double *upper, double *avnode, 
				double *bestcrit, int *nodestatus, 
				int *nodepop, int *nodestart, int *nrnodes, 
				int *nthsize, double *rsnodecost, int *ncase, 
				int *parent, double *ut, double *v, 
				double *xt, int *mtry, int *ip, int *mbest, 
				int *cat, double *tgini, int *mind);
extern void F77_NAME(rtestreebag)(double *x, int *nsample, int *mdim, 
				 int *treemap, int *nodestatus, int *nrnodes,
				 int *ndbigtree, double *ytr, double *upper, 
				 double *avnode, int *mbest, int *cat);
