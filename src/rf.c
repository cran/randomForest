/*****************************************************************
   Copyright (C) 2001-2 Leo Breiman, Adele Cutler, Andy Liaw and Matthew Wiener
  
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

   C driver for Breiman & Cutler's random forest code.
   Re-written from the original main program in Fortran.
   Andy Liaw Feb. 7, 2002.
   Modifications to get the forest out Matt Wiener Feb. 26, 2002. 
 *****************************************************************/

#include <stdio.h>
#include <assert.h>
#include <R.h>
#include "rf.h"

/*  Define the R RNG for use from Fortran. */
void F77_SUB(rrand)(double *r) { *r = unif_rand(); }

void rf(double *x, int *ncol, int *nrow, int *cl, int *ncl, int *cat, 
        int *maxcat, int *addcl, int *ntree, int *nvar, int *ipi, double
        *pi, int *nodesize, int *imp, int *iprox, int *noutlier, 
        double *outlier, int *outcl, int *counttr, double *prox, 
	double *imprt, int *trace, int *ndbigtree, int *nodestatus, 
	int *bestvar, int *treemap, int *nodeclass, double *xbestsplit, 
	double *pid, int *keepf, int *testdat, double *xts, int *clts, int *nts, 
	double *countts, int *outclts, int *labelts)
{
  /******************************************************************
   *  C wrapper for random forests:  get input from R and drive
   *  the Fortran routines.
   *
   *  Input:
   *
   *  x:        matrix of predictors (transposed!)
   *  ncol:     number of predictor variables
   *  nrow:     number of cases in the training set
   *  cl:       class labels of the data
   *  ncl:      number of classes in the response
   *  cat:      integer vector of number of classes in the predictor;
   *            1=continuous
   *  addcl:    =0 if the data has class labels;
   *            =1 add a synthetic class by sampling from product 
   *               of marginals
   *            =2 add a synthetic class by sampling from product
   *               of uniforms
   *  ntree:    number of trees
   *  nvar:     number of predictors to use for each split
   *  ipi:      0=use class proportion as prob.; 1=use supplied priors
   *  pi:       double vector of class priors
   *  nodesize: minimum node size: no node with fewer than ndsize
   *            cases will be split
   *  imp:      =0 turns off variable importance; =1 turns it on
   *  iprox:    =0 turns off proximity calculation; =1 turns it on
   *  noutlier: =0 turns off outlyingness measures; =1 turns it on
   *
   *  Output:
   *
   *  outcl:    class predicted by RF
   *  counttr:  matrix of votes (transposed!)
   *  imprt:    matrix of variable importance measures (if imp=1)
   *  prox:     matrix of proximity (if iprox=1)
   *  outlier:  measure of outlyingness (if noutlier=1)
   ******************************************************************/

  int nsample0, mdim, nclass, iaddcl, jbt, mtry, ntest, nsample, ndsize,
    nrnodes, mimp, nimp, near, nuse;
  int jb, j, n, m, mr, k, kt, itwo, arrayindex;

  int *out, *bestsplitnext, *bestsplit,
    *nodepop, *parent, *jin, *ndble, *nodex,
    *nodexts, *nodestart, *ta, *ncase, *jerr, *iv, *isort, *ncp, *clp,
    *jtr, *nc, *msum, *idmove, *jvr, *countimp,
    *at, *a, *b, *cbestsplit, *mind, *jts;
  
  double errtr, errts, errc;

  double *tgini, *v, *tx, *wl, *classpop, *errimp,
    *rimpmarg, *tclasscat, *tclasspop, *rmargin, *win, *tp,
    *wc, *wr, *wtt, *diffmarg, *cntmarg, *rmissimp,
    *tout, *tdx, *sm, *p, *q, *iw;

  nsample0 = *nrow;
  mdim     = *ncol;
  nclass   = (*ncl==1) ? 2 : *ncl;
  iaddcl   = *addcl;
  ndsize   = *nodesize;
  jbt      = *ntree;
  mtry     = *nvar;
  ntest    = *nts; 
  nsample = (iaddcl > 0) ? (nsample0 + nsample0) : nsample0;
  nrnodes = 2 * (nsample / ndsize) + 1;
  mimp = (*imp == 1) ? mdim : 1;
  nimp = (*imp == 1) ? nsample : 1;
  near = (*iprox == 1) ? nsample0 : 1;
  if (*trace == 0) *trace = jbt + 1;

  tgini =      (double *) S_alloc(mdim, sizeof(double));
  v =          (double *) S_alloc(nsample, sizeof(double));
  tx =         (double *) S_alloc(nsample, sizeof(double));
  wl =         (double *) S_alloc(nclass, sizeof(double));
  classpop =   (double *) S_alloc(nclass*nrnodes, sizeof(double));
  errimp =     (double *) S_alloc(mimp, sizeof(double));
  rimpmarg =   (double *) S_alloc(mdim*nsample, sizeof(double));
  tclasscat =  (double *) S_alloc(nclass*32, sizeof(double));
  tclasspop =  (double *) S_alloc(nclass, sizeof(double));
  rmargin =    (double *) S_alloc(nsample, sizeof(double));
  win =        (double *) S_alloc(nsample, sizeof(double));
  tp =         (double *) S_alloc(nsample, sizeof(double));
  wc =         (double *) S_alloc(nclass, sizeof(double));
  wr =         (double *) S_alloc(nclass, sizeof(double));
  wtt =        (double *) S_alloc(nsample, sizeof(double));
  diffmarg =   (double *) S_alloc(mdim, sizeof(double));
  cntmarg =    (double *) S_alloc(mdim, sizeof(double));
  rmissimp =   (double *) S_alloc(mimp, sizeof(double));
  tout =       (double *) S_alloc(near, sizeof(double));
  tdx =        (double *) S_alloc(nsample0, sizeof(double));
  sm =         (double *) S_alloc(nsample0, sizeof(double));
  p =          (double *) S_alloc(nsample0, sizeof(double));
  q =          (double *) S_alloc(nclass*nsample, sizeof(double));
  iw =         (double *) S_alloc(nsample, sizeof(double));
  
  out =           (int *) S_alloc(nsample, sizeof(int));
  countimp =      (int *) S_alloc(nclass*nimp*mimp, sizeof(int));
  bestsplitnext = (int *) S_alloc(nrnodes, sizeof(int));
  bestsplit =     (int *) S_alloc(nrnodes, sizeof(int));
  nodepop =       (int *) S_alloc(nrnodes, sizeof(int));
  parent =        (int *) S_alloc(nrnodes, sizeof(int));
  jin =           (int *) S_alloc(nsample, sizeof(int));
  ndble =         (int *) S_alloc(nsample0, sizeof(int));
  nodex =         (int *) S_alloc(nsample, sizeof(int));
  nodexts =       (int *) S_alloc(ntest, sizeof(int));
  nodestart =     (int *) S_alloc(nrnodes, sizeof(int));
  ta =            (int *) S_alloc(nsample, sizeof(int));
  ncase =         (int *) S_alloc(nsample, sizeof(int));
  jerr =          (int *) S_alloc(nsample, sizeof(int));
  iv =            (int *) S_alloc(mdim, sizeof(int)); 
  isort =         (int *) S_alloc(nsample, sizeof(int));
  ncp =           (int *) S_alloc(near, sizeof(int));
  clp =           (int *) S_alloc(near, sizeof(int));
  jtr =           (int *) S_alloc(nsample, sizeof(int));
  nc =            (int *) S_alloc(nclass, sizeof(int));
  msum =          (int *) S_alloc(mdim, sizeof(int));
  jts =           (int *) S_alloc(ntest, sizeof(int));
  idmove =        (int *) S_alloc(nsample, sizeof(int));
  jvr =           (int *) S_alloc(nsample, sizeof(int));
  at =            (int *) S_alloc(mdim*nsample, sizeof(int));
  a =             (int *) S_alloc(mdim*nsample, sizeof(int));
  b =             (int *) S_alloc(mdim*nsample, sizeof(int));
  cbestsplit =    (int *) S_alloc(*maxcat*nrnodes, sizeof(int));
  mind =          (int *) S_alloc(mdim, sizeof(int));
  
  /* SET UP DATA TO ADD A CLASS++++++++++++++++++++++++++++++ */

  if(iaddcl >= 1) 
    F77_CALL(createclass)(x, cl, &nsample0, &nsample, &mdim, tdx, p,
			  sm, ndble, &iaddcl);
  
  
  /*    INITIALIZE FOR RUN */
  if(*testdat == 1) {
    for (n = 0; n < nclass; ++n) {
      for (k = 0; k < ntest; ++k) countts[n + k*nclass] = 0.0;
    }
  }

  F77_CALL(zerm)(counttr, &nclass, &nsample);
  F77_CALL(zerv)(out, &nsample);
  F77_CALL(zervr)(tgini, &mdim);
  F77_CALL(zerv)(msum, &mdim);

  errtr = 0.0;  errts = 0.0;
  
  if(*imp == 1) {
    for(m = 0; m < mdim; ++m) {
      for(k = 0; k < nsample; ++k) {
	for(j = 0; j < nclass; ++j) {
	  assert(j + (k*nclass) + (m*nclass*nsample) < nclass*nimp*mimp);
	  countimp[j + (k*nclass) + (m*nclass*nsample)] = 0;
	}
      }
    }
  }
  
  F77_CALL(prep)(cl, &nsample, &nclass, ipi, pi, pid, nc, wtt);
  F77_CALL(makea)(x, &mdim, &nsample, cat, isort, v, at, b, &mdim);


  /*   START RUN   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  itwo = 2;
  GetRNGstate(); 

  for(jb = 0; jb < jbt; jb++) {

    arrayindex = (*keepf == 1) ? jb * nrnodes : 0;

    F77_CALL(zerv)(nodestatus + arrayindex, &nrnodes);
    F77_CALL(zerm)(treemap + 2*arrayindex, &itwo, &nrnodes);
    F77_CALL(zervr)(xbestsplit + arrayindex, &nrnodes);
    F77_CALL(zerv)(nodeclass + arrayindex, &nrnodes);
    F77_CALL(zerv)(jin, &nsample);
    F77_CALL(zervr)(tclasspop, &nclass);
    F77_CALL(zervr)(win, &nsample);

    for(n = 0; n < nsample; n++) {
      k = (unif_rand() * nsample) - 1;
      assert(k < nsample); 
      assert(k >= 0); 
      tclasspop[cl[k] - 1] = tclasspop[cl[k] - 1] + wtt[k];
      win[k] = win[k] + wtt[k];
      jin[k] = 1;
    }
    
    F77_CALL(eqm)(a, at, &mdim, &nsample);
    F77_CALL(moda)(a, &nuse, &nsample, &mdim, cat, maxcat, ncase, jin, ta);
    
    F77_CALL(buildtree)(a, b, cl, cat, &mdim, &nsample, &nclass, 
			treemap + 2*arrayindex, 
			bestvar + arrayindex, bestsplit, bestsplitnext, tgini,
			nodestatus + arrayindex, nodepop, nodestart, classpop,
			tclasspop, tclasscat, ta, &nrnodes, idmove,
			&ndsize, ncase, parent, jin, &mtry, iv,
			nodeclass + arrayindex, ndbigtree + jb, win, wr, wc,
			wl, &mdim, &nuse, mind); 

    F77_CALL(xtranslate)(x, &mdim, &nrnodes, &nsample, bestvar + arrayindex, 
			 bestsplit, bestsplitnext, xbestsplit + arrayindex,
			 nodestatus + arrayindex, cat, ndbigtree + jb); 

    /*  Get test set error */
    if(*testdat == 1) {
      F77_CALL(testreebag)(xts, &ntest, &mdim, treemap, nodestatus, 
			   xbestsplit, cbestsplit, bestvar, nodeclass,
			   &nrnodes, ndbigtree + jb, cat, &nclass,
			   jts,nodexts,maxcat);
      F77_CALL(comptserr)(countts, jts, clts, outclts, &ntest, &nclass,
			  &errts, pid, labelts);
    }
    
    /*  GET OUT-OF-BAG ESTIMATES */
    F77_CALL(testreebag)(x, &nsample, &mdim, treemap + 2*arrayindex, 
			 nodestatus + arrayindex, xbestsplit + arrayindex, 
			 cbestsplit, bestvar + arrayindex, 
			 nodeclass + arrayindex, &nrnodes, ndbigtree + jb, 
			 cat, &nclass, jtr, nodex, maxcat); 

    for(n = 0; n < nsample; n++) {
      if(jin[n] == 0) {
	assert((n*nclass + jtr[n] - 1) < (*ncl * nsample)); 
	assert((n*nclass + jtr[n] - 1) >= 0); 
	counttr[n*nclass + jtr[n] - 1] ++;
        out[n]++;
      }
    }

    if ((jb + 1) % *trace == 0 || jb + 1 == jbt) {
      F77_CALL(oob)(&nsample, &nclass, jin, cl, jtr, jerr, counttr, out,
		    &errtr, &errc, rmargin, q, outcl, wtt);
      if(*trace < jbt) {
	Rprintf("%4i: OOB error rate=%5.2f%%\t", jb+1, 100.0*errtr);
	if(*labelts == 1) 
 	  Rprintf("Test set error rate=%5.2f%%", 100.0*errts);
 	Rprintf("\n");
      }
    }

    /*  DO VARIABLE IMPORTANCE  */
    if(*imp == 1) {
      F77_CALL(zerv)(iv, &mdim);
      for(kt = 0; kt < *(ndbigtree + jb); kt++) {
        if(nodestatus[kt + arrayindex] != -1) {
	  assert(kt + arrayindex < jbt * nrnodes);
	  assert(kt + arrayindex >= 0);
	  assert(bestvar[kt + arrayindex] - 1 < mdim);
	  assert(bestvar[kt + arrayindex] - 1 >= 0);
	  iv[bestvar[kt + arrayindex] - 1] = 1;
	  msum[bestvar[kt + arrayindex] - 1] ++;
	}
      }
      for(mr = 1; mr <= mdim; mr++) {
        if(iv[mr-1] == 1) {
	  F77_CALL(permobmr)(&mr, x, tp, tx, jin, &nsample, &mdim);
	  F77_CALL(testreebag)(x, &nsample, &mdim, treemap + 2*arrayindex, 
			       nodestatus + arrayindex,
			       xbestsplit + arrayindex, cbestsplit, 
			       bestvar + arrayindex, nodeclass + arrayindex, 
			       &nrnodes, ndbigtree + jb, cat,
			       &nclass, jvr, nodex, maxcat); 
	  for(n = 0; n < nsample; n++) {
	    if(jin[n] == 0 && jtr[n] != jvr[n]) {
	      assert(jvr[n]-1 + n*nclass + (mr-1)*nclass*nsample < 
		     nclass * nimp * mimp);
	      assert(jvr[n]-1 + n*nclass + (mr-1)*nclass*nsample >= 0);
	      countimp[jvr[n]-1 + n*nclass + (mr-1)*nclass*nsample]++;
	      countimp[jtr[n]-1 + n*nclass + (mr-1)*nclass*nsample]--;
	    }
	  }
	  for(n = 0; n < nsample; n++) {
	    assert(mr-1 + n*mdim < mdim*nsample);
	    assert(mr-1 + n*mdim >= 0);
	    x[mr-1 + n*mdim] = tx[n];
	  }
	}
      }
    }
        
    /*  DO PROXIMITIES */
    if(*iprox == 1) {
      for(n = 0; n < near; n++) {
	for(k = 0; k < near; k++) {
	  assert(k < nsample);
	  assert(n < nsample);
	  assert(k >= 0);
	  assert(n >= 0);
	  assert(k*near + n < near * near);
	  assert(k*near + n >= 0);
	  if(nodex[k] == nodex[n]) prox[k*near + n] += 1.0;
	}
      }
    }
  }
  PutRNGstate();
      
  /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/  
  if(*iprox == 1) {
    for(n = 0; n < near; n++) {
      for(k = n + 1; k < near; k++) {
	assert(k*near + n < near * near);
	assert(n*near + k < near * near);
	assert(k*near + n >= 0);
	assert(n*near + k >= 0);
	prox[near*k + n] = prox[near*k + n] / jbt;
	prox[near*n + k] = prox[near*k + n];
      }
      assert(near*n + n < near * near);
      assert(near*n + n >= 0);
      prox[near*n + n] = 1.0;
    }
  }

  /*  OUTLIER DATA +++++++++++++++++++++++++++++*/
  if(*noutlier == 1) {
    F77_CALL(locateout)(prox, cl, &near, &nsample, &nclass, ncp, &iaddcl,
			outlier, tout, isort, clp);
  }
  
  /*  IMP DATA ++++++++++++++++++++++++++++++++++*/
  if (*imp == 1) {
    F77_CALL(finishimp)(rmissimp, countimp, out, cl, &nclass, &mdim,
			&nsample, errimp, rimpmarg, diffmarg, cntmarg, 
			rmargin, counttr, outcl, &errtr);

  /*	GIVES STANDARD IMP OUTPUT */
    for(m = 0; m < mdim; m++) {
      assert(m < mimp);
      assert(m >= 0);
      assert(3*mdim + m < 4*mdim);
      assert(3*mdim + m >= 0);
      imprt[m] = errimp[m];
      imprt[mdim + m] = ((diffmarg[m] > 0.0) ? diffmarg[m] : 0.0);
      imprt[2*mdim + m] = ((cntmarg[m] > 0.0) ? cntmarg[m] : 0.0);
      imprt[3*mdim + m] = tgini[m] / jbt;
    }
  }
}
