/*******************************************************************
   Copyright (C) 2001-4 Leo Breiman, Adele Cutler and Merck & Co., Inc.
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.                            
*******************************************************************/

void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat, 
	int *sampsize, int *Options, int *ntree, int *nvar,
	int *ipi, double *pi, double *cut, int *nodesize, 
        int *outcl, int *counttr, double *prox, 
	double *imprt, double *, double *impmat, int *nrnodes, int *ndbigtree, 
	int *nodestatus, int *bestvar, int *treemap, int *nodeclass,
	double *xbestsplit, double *pid, double *errtr, 
	int *testdat, double *xts, int *clts, int *nts, double *countts,
	int *outclts, int *labelts, double *proxts, double *errts);

void runforest(int *mdim, int *ntest, int *nclass, int *maxcat,
	       int *nrnodes, int *jbt,
	       double *xts, double *xbestsplit, double *pid, 
	       double *cutoff, double *countts, int *treemap,
	       int *nodestatus, int *cat, int *cbestsplit,
	       int *nodeclass, int *jts, int *jet, int *bestvar,
	       int *nodexts, int *ndbigtree, int *keepPred, 
	       int *prox, double *proxmatrix, int *nodes);

void regTree(double *x, double *y, int mdim, int nsample, 
	     int *treemap, double *upper, double *avnode, int *nodestatus, 
	     int nrnodes, int nthsize, int mtry, int *mbest, int *cat, 
	     double *tgini, int *varUsed);

void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample, 
		   int ndstart, int ndend, int *msplit, double *decsplit, 
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat);

void predictRegTree(double *x, int nsample, int mdim, int *doPred,
		    int *treemap, int *nodestatus, int nrnodes, 
		    int ndbigtree, double *ypred, double *split, 
		    double *nodepred, int *bestvar, int *cat, int *nodex);

void predictClassTree(double *x, int n, int mdim, int *doPred, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *cbestsplit, int *bestvar, int *nodeclass,
		      int nrnodes, int ndbigtree, int *cat, int nclass,
		      int *jts, int *nodex, int maxcat);

double pack(int l, int *icat);
void unpack(int l, int npack, int *icat);

void zeroInt(int *x, int length);
void zeroDouble(double *x, int length);
void createClass(double *x, int realN, int totalN, int mdim);
void prepare(int *cl, const int nsample, const int nclass, const int ipi, 
	     double *pi, double *pid, int *nc, double *wtt);
void makeA(double *x, const int mdim, const int nsample, int *cat, int *a, 
           int *b);
void modA(int *a, int *nuse, const int nsample, const int mdim, int *cat, 
          const int maxcat, int *ncase, int *jin);
void Xtranslate(double *x, int mdim, int nrnodes, int nsample, 
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize);
void permuteOOB(int m, double *x, int *in, int nsample, int mdim);

/* Template of Fortran subroutines to be called from the C wrapper */
extern void F77_NAME(zerm)(int *,     int *, int *);
extern void F77_NAME(zermd)(double *, int *, int *);
extern void F77_NAME(zerv)(int *,     int *);
extern void F77_NAME(zervr)(double *, int *);
extern void F77_NAME(buildtree)(int *, int *, int *, int *, int *, int *, int
				*, int *, int *, int *, int *, double *, int
				*, int *, int *, double *, double *, double
				*, int *, int *, int *, int *, int *, int *,
				int *, int *, int *, int *, int *, double *,
				double *, double *, double *, int *, int *,
				int *); 
