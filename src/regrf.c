/************************************************************************
     Last change:  LB   13 Mar 2002   12:01 pm

      copyright 1999 by leo Breiman
      this is free software and can be used for any purpose. 
      It comes with no guarantee.  
 *************************************************************************/
#include <R.h>
#include <Rmath.h>
#include "rf.h"

/*  Define the R RNG for use from Fortran. */
/* void F77_SUB(rrand)(double *r) { *r = unif_rand(); } */

void regrf(double *x, double *y, int *nsample, int *mdim, int *nthsize, 
	   int *nrnodes, int *jbt, int *mtry, int *imp, int *cat, int *jprint,
	   int *iprox, double *yptr, double *errimp, double *prox, 
	   int *ndbigtree, int *nodestatus, int *treemap, double *avnode, 
	   int *mbest, double *upper, double *mse, double *rsq, int *keepf, 
	   int *testdat, double *xts, int *nts, double *yts, int *labelts,
	   double *ypred, double *proxts, double *msets)
{
  /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases

   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   jbt=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be 
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/
  
  double errts = 0.0, averrb, avy, avyts, vary, varyts, astr, asd, xrand, 
    em, errb = 0.0;

  double *yb, *rsnodecost, *bestcrit, *sd, *wts, *v, *ut, *xt, *xb, *ytr, 
    *yl, *xa, *utr, *predimp, *za, *ytree, *tgini, maximp;
  
  int i, k, m, mr, mrind, n, nls, ntrue, jout, nimp, mimp, jb, idx, ntest;
  
  int *jdex, *nodepop, *npert, *ip, *nperm, *parent, *nout, *jin, *isort, 
    *nodestart, *ncase, *nbrterm, *jperm, *incl, *mind, *nodex, *nodexts;
  
  nimp = (*imp == 1) ? (*imp * *nsample) : 1;
  mimp = (*imp == 1) ? (*imp * *mdim) : 1;
  ntest = *nts;

  if(*jprint == 0) *jprint = *jbt + 1;

  yb         = (double *) R_alloc(*nsample, sizeof(double));
  rsnodecost = (double *) R_alloc(*nrnodes, sizeof(double));
  bestcrit   = (double *) R_alloc(*nrnodes, sizeof(double));
  sd         = (double *) R_alloc(*mdim, sizeof(double));
  wts        = (double *) R_alloc(*nsample, sizeof(double));
  v          = (double *) R_alloc(*nsample, sizeof(double));
  ut         = (double *) R_alloc(*nsample, sizeof(double));
  xt         = (double *) R_alloc(*nsample, sizeof(double));
  xb         = (double *) R_alloc(*mdim * *nsample, sizeof(double));
  ytr        = (double *) R_alloc(*nsample, sizeof(double));
  ytree      = (double *) R_alloc(ntest, sizeof(double));
  yl         = (double *) R_alloc(*nsample, sizeof(double));
  xa         = (double *) R_alloc(3 * *mdim, sizeof(double));
  utr        = (double *) R_alloc(*nsample, sizeof(double));
  predimp    = (double *) R_alloc(nimp * mimp, sizeof(double));
  za         = (double *) R_alloc(*mdim, sizeof(double));
  tgini      = (double *) R_alloc(*mdim, sizeof(double));

  jdex       = (int *) R_alloc(*nsample, sizeof(int));
  nodepop    = (int *) R_alloc(*nrnodes, sizeof(int));
  npert      = (int *) R_alloc(*nsample, sizeof(int));
  ip         = (int *) R_alloc(*mdim, sizeof(int));
  nperm      = (int *) R_alloc(*nsample, sizeof(int));
  parent     = (int *) R_alloc(*nrnodes, sizeof(int));
  nout       = (int *) R_alloc(*nsample, sizeof(int));
  jin        = (int *) R_alloc(*nsample, sizeof(int));
  isort      = (int *) R_alloc(*nsample, sizeof(int));
  nodestart  = (int *) R_alloc(*nrnodes, sizeof(int));
  ncase      = (int *) R_alloc(*nsample, sizeof(int));
  nbrterm    = (int *) R_alloc(*nrnodes, sizeof(int));
  jperm      = (int *) R_alloc(*jbt, sizeof(int));
  incl       = (int *) R_alloc(*mdim, sizeof(int));
  mind       = (int *) R_alloc(*mdim, sizeof(int)); 
  nodex      = (int *) R_alloc(*nsample, sizeof(double));
  if(*testdat == 1)  nodexts = (int *) R_alloc(ntest, sizeof(double));

    
  averrb = 0.0;
	
  avy = 0.0;
  vary = 0.0;

  for(n=0; n < *nsample; ++n) {
    yptr[n] = 0.0;
    nout[n] = 0;
    ntrue = n;
    vary += ntrue * (y[n] - avy)*(y[n] - avy) / (ntrue + 1);
    avy = (ntrue * avy + y[n]) / (ntrue + 1);
  }
  vary /= *nsample;

  varyts = 0.0;
  avyts = 0.0;
  if(*testdat == 1) {
    for(n = 0; n < ntest; ++n) {
      ntrue = n;
      varyts += ntrue * (yts[n] - avyts)*(yts[n] - avyts) / (ntrue + 1);
      avyts = (ntrue * avyts + yts[n]) / (ntrue + 1);
    }
    varyts /= ntest;
  }

  astr = 0.0;
  asd = 0.0;

  if(*iprox == 1) {
    for(n = 0; n < *nsample; ++n) {
      for(k = 0; k < *nsample; ++k) 
	prox[k * *nsample + n] = 0.0;
    }
    if(*testdat == 1) {
      for(n = 0; n < ntest; ++n) {
	for(k = 0; k < ntest + *nsample; ++k) 
	  prox[k * ntest + n] = 0.0;
      }
    }
  }

  for(m = 0; m < *mdim; ++m) {
    tgini[m] = 0.0;
  }
  
  if(*imp==1) {
    for(m = 0; m < *mdim; ++m) {
      errimp[m] = 0.0;
      errimp[m + *mdim] = 0.0;
      for(n = 0; n < *nsample; ++n) {
	predimp[n + m * *nsample] = 0.0;
      }
    }
  }

  if(*labelts == 1) {
    for(n = 0; n < ntest; ++n) ypred[n] = 0.0;
  }

  GetRNGstate();

  /*************************************
   Start the loop over trees.
  *************************************/
  for(jb = 0; jb < *jbt; ++jb) {
    idx = (*keepf == 1) ? jb * *nrnodes : 0;
    for(n = 0; n < *nsample; ++n) {
      jin[n] = 0;
      nodex[n] = 0;
    }

    for(n = 0; n < *nsample; ++n) {
      xrand = unif_rand();
      k = xrand * *nsample;
      jin[k] = 1;
      yb[n] = y[k];
      for(m = 0; m < *mdim; ++m) {
	xb[m + n * *mdim] = x[m + k * *mdim];
      }
    }

    nls = *nsample;
    
    F77_CALL(rbuildtree)(xb, yb, yl, mdim, &nls, nsample, treemap + (2*idx), 
			jdex, upper + idx, avnode + idx, bestcrit, 
			nodestatus + idx, nodepop,
			nodestart, nrnodes, nthsize, rsnodecost, ncase, 
			parent, ut, v, xt, mtry, ip, mbest + idx, cat, 
			tgini, mind);

    ndbigtree[jb] = *nrnodes;
    for(k = *nrnodes-1; k >= 0; --k) {
      if (nodestatus[k + idx]==0) ndbigtree[jb]--;
      if (nodestatus[k + idx]==2) nodestatus[k + idx] = -1;
    }

    for(n = 0; n < *nsample; ++n) ytr[n] = 0.0;
    
    F77_CALL(rtestreebag)(x, nsample, mdim, treemap + 2*idx, nodestatus + idx, 
			 nrnodes, ndbigtree + jb, ytr, upper + idx, 
			 avnode + idx, mbest + idx, cat, nodex);
    
    errb = 0.0;
    jout = 0;
    for (n = 0; n < *nsample; ++n) {
      if (jin[n] == 0) { 
	yptr[n] = (nout[n] * yptr[n] + ytr[n]) / (nout[n] + 1);
	nout[n] ++;
      }
      
      if (nout[n] > 0) jout ++;
      errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
    }
    errb /= *nsample;

    if(*testdat == 1) {
      for(i = 0; i < ntest; ++i) {
	ytree[i] = 0.0;
	nodexts[i] = 0;
      }
      F77_CALL(rtestreebag)(xts, &ntest, mdim, treemap + 2*idx, 
			    nodestatus + idx, nrnodes, ndbigtree + jb,
			    ytree, upper + idx, avnode + idx, 
			    mbest + idx, cat, nodexts);
      errts = 0.0;
      for(n = 0; n < ntest; ++n) {
	ypred[n] = (jb * ypred[n] + ytree[n]) / (jb + 1);
	if(*labelts == 1)
	  errts += (yts[n] - ypred[n]) * (yts[n] - ypred[n]);
      }
      if(*labelts == 1) errts /= ntest;
    }

    if ((jb + 1) % *jprint == 0) {
      Rprintf("%d: ", jb + 1);
      if(*labelts == 1) Rprintf("MSE(Test)=%f  %%Var(y)=%7.2f  ", 
				errts, 100.0 * errts / varyts);
      Rprintf("MSE(OOB)=%f  %%Var(y)=%7.2f\n", errb, 100*errb/vary);
    }
    mse[jb] = errb;
    if(*labelts == 1) msets[jb] = errts;

    /*  DO PROXIMITIES */
    if(*iprox == 1) {
      for(n = 0; n < *nsample; ++n) {
	for(k = 0; k < *nsample; ++k) {
	  if(nodex[k] == nodex[n]) prox[k * *nsample + n] += 1.0;
	}
      }

      /* proximity for test data */
      if(*testdat == 1) {
	for(n = 0; n < ntest; ++n) {
	  for(k = 0; k <= n; ++k) {
	    if(nodexts[k] == nodexts[n]) {
	      proxts[k * ntest + n] += 1.0;
	      proxts[n * ntest + k] = proxts[k * ntest + n];
	    }
	  }
	  for(k = 0; k < *nsample; ++k) {
	    if(nodexts[n] == nodex[k]) proxts[n + ntest * (k+ntest)] += 1.0; 
	  } 
	}
      } 
    }

    /* Variable importance */
    if(*imp == 1) { 
      for (mr = 0; mr < *mdim; ++mr) {
	mrind = mr + 1;
	F77_CALL(permobmr)(&mrind, x, utr, xt, jin, nsample, mdim);
	F77_CALL(rtestreebag)(x, nsample, mdim, treemap + 2*idx, 
			     nodestatus + idx, nrnodes, 
			     ndbigtree + jb, ytr, upper + idx, avnode + idx, 
			     mbest + idx, cat, nodex);
	for (n = 0; n < *nsample; ++n) {
	  x[mr + n * *mdim] = xt[n];
	}
	em = 0.0;
	for (n = 0; n < *nsample; ++n) {
	  if(jin[n] == 0) { 
	    predimp[n + mr * *nsample] = (nout[n] * predimp[n + mr * *nsample]
					  + ytr[n]) / (nout[n] + 1);
	  }
	  em += (y[n] - predimp[n + mr* *nsample]) *
	    (y[n] - predimp[n + mr* *nsample]);
	  /*	  if(mr==0) Rprintf("%i %10.5f\n", n+1, predimp[n]);*/
	}
	errimp[mr] = em / *nsample;
      }
    }
  }
  PutRNGstate();

  /* end of tree iterations=======================================*/
  *rsq = 1.0 - errb/vary;
  maximp = 0.0;
  for(m = 0; m < *mdim; ++m) {
    if(tgini[m] > maximp) maximp = tgini[m];
  }

  if (*imp == 1) {
    for (m = 0; m < *mdim; ++m) {
      errimp[m] = 100 * ((errimp[m] / errb) - 1);
      if(errimp[m] <= 0.0) errimp[m] = 0.0;
      errimp[m + *mdim] = tgini[m] / maximp;
    }
  } else {
    for(m = 0; m < *mdim; ++m) 
      errimp[m] = tgini[m] / maximp;
  }
  
}


void runrforest(double *xts, double *ypred, int *mdim, int *ntest, int *ntree, 
		int *ndbigtree, int *treemap, int *nodestatus, int *nrnodes, 
		double *upper, double *avnodes, int *mbest, int *cat,
		int *iprox, double *proximity) 
{
  
  int i, j, k, n, idx, *nodex;
  double *ytree;

  ytree = (double *) R_alloc(*ntest, sizeof(double));
  nodex = (int *) R_alloc(*ntest, sizeof(int));

  for(i = 0; i < *ntest; ++i) {
    ytree[i] = 0.0;
    if(*iprox == 1) 
      for(j = 0; j < *ntest; ++j) {
	proximity[i + j * *ntest] = 0.0;
	proximity[j + i * *ntest] = 0.0;
      }
  }

  for(i = 0; i < *ntree; ++i) {
    idx = i * *nrnodes;
    for(j = 0; j < *ntest; j++) ytree[j] = 0.0;

    F77_CALL(rtestreebag)(xts, ntest, mdim, treemap + 2*idx, nodestatus + idx,
			 nrnodes, ndbigtree + i, ytree, upper + idx, 
			 avnodes + idx, mbest + idx, cat, nodex);

    for(j = 0; j < *ntest; ++j) ypred[j] += ytree[j];

    /* if desired, do proximities for this round */
    if (*iprox == 1) {
      for(n = 0; n < *ntest; ++n) {
	for(k = 0; k <= n; ++k) {
	  if(nodex[n] == nodex[k]) {
	    proximity[n + *ntest * k] += 1.0;
	    proximity[k + *ntest * n] = proximity[n + *ntest * k];
	  }
	}
      }
    }
  }
  for(i = 0; i < *ntest; ++i) ypred[i] /= *ntree;

  if(*iprox == 1) {
    for(i = 0; i < *ntest; ++i) {
      for(j = 0; j < *ntest; ++j) {
	proximity[i * *ntest + j] /= *ntree;
      }
    }
  }
}
