/************************************************************************
      copyright 1999 by leo Breiman
      this is free software and can be used for any purpose. 
      It comes with no guarantee.  
 *************************************************************************/
#include <R.h>
#include "rf.h"

void simpleLinReg(int nsample, double *x, double *y, double *coef, 
		  double *mse, int *hasPred);


void regrf(double *x, double *y, int *nsample, int *mdim, int *nthsize,
	   int *nrnodes, int *jbt, int *mtry, int *imp, int *cat, 
	   int *jprint, int *iprox, int *oobprox, int *biasCorr, double *yptr, 
	   double *errimp, double *prox, int *ndbigtree, int *nodestatus,
	   int *treemap, double *avnode, int *mbest, double *upper,
	   double *mse, double *rsq, int *keepf, int *testdat,
	   double *xts, int *nts, double *yts, int *labelts, 
	   double *ypred, double *proxts, double *msets, double *coef, 
	   int *nout)
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
    errb = 0.0, resid=0.0, ooberr, ooberrperm, delta;

  double *yb, *rsnodecost, *bestcrit, *sd, *wts, *v, *ut, *xt, *xb, *ytr, 
    *yl, *xa, *utr, *za, *ytree, *tgini, *sse; /* *predimp, */
  
  int i, k, m, mr, mrind, n, nls, ntrue, jout, nimp, mimp, jb, idx, ntest;
  int *oobpair;

  int *jdex, *nodepop, *npert, *ip, *nperm, *parent, *jin, *isort, 
    *nodestart, *ncase, *nbrterm, *jperm, *incl, *mind, *nodex, *nodexts;
  
  nimp = (*imp) ? *nsample : 1;
  mimp = (*imp) ? *mdim : 1;
  ntest = *nts;

  if (*jprint == 0) *jprint = *jbt + 1;

  yb         = (double *) S_alloc(*nsample, sizeof(double));
  wts        = (double *) S_alloc(*nsample, sizeof(double));
  v          = (double *) S_alloc(*nsample, sizeof(double));
  ut         = (double *) S_alloc(*nsample, sizeof(double));
  xt         = (double *) S_alloc(*nsample, sizeof(double));
  ytr        = (double *) S_alloc(*nsample, sizeof(double));
  yl         = (double *) S_alloc(*nsample, sizeof(double));
  utr        = (double *) S_alloc(*nsample, sizeof(double));
  rsnodecost = (double *) S_alloc(*nrnodes, sizeof(double));
  bestcrit   = (double *) S_alloc(*nrnodes, sizeof(double));
  xa         = (double *) S_alloc(3 * *mdim, sizeof(double));
  sd         = (double *) S_alloc(*mdim, sizeof(double));
  xb         = (double *) S_alloc(*mdim * *nsample, sizeof(double));
  sse        = (double *) S_alloc(mimp, sizeof(double));
  za         = (double *) S_alloc(*mdim, sizeof(double));
  tgini      = (double *) S_alloc(*mdim, sizeof(double));

  jdex       = (int *) S_alloc(*nsample, sizeof(int));
  npert      = (int *) S_alloc(*nsample, sizeof(int));
  nperm      = (int *) S_alloc(*nsample, sizeof(int));
  jin        = (int *) S_alloc(*nsample, sizeof(int));
  isort      = (int *) S_alloc(*nsample, sizeof(int));
  ncase      = (int *) S_alloc(*nsample, sizeof(int));
  nodex      = (int *) S_alloc(*nsample, sizeof(int));
  nodepop    = (int *) S_alloc(*nrnodes, sizeof(int));
  parent     = (int *) S_alloc(*nrnodes, sizeof(int));
  nodestart  = (int *) S_alloc(*nrnodes, sizeof(int));
  nbrterm    = (int *) S_alloc(*nrnodes, sizeof(int));
  incl       = (int *) S_alloc(*mdim, sizeof(int));
  mind       = (int *) S_alloc(*mdim, sizeof(int)); 
  ip         = (int *) S_alloc(*mdim, sizeof(int));
  jperm      = (int *) S_alloc(*jbt, sizeof(int));
  if (*testdat) {
    ytree      = (double *) S_alloc(ntest, sizeof(double));
    nodexts    = (int *) S_alloc(ntest, sizeof(int));
  }
  if (*oobprox) oobpair = (int *) S_alloc(*nsample * *nsample, sizeof(int));

  averrb = 0.0;
  avy = 0.0;
  vary = 0.0;

  for (n = 0; n < *nsample; ++n) {
    yptr[n] = 0.0;
    nout[n] = 0;
    ntrue = n;
    vary += ntrue * (y[n] - avy)*(y[n] - avy) / (ntrue + 1);
    avy = (ntrue * avy + y[n]) / (ntrue + 1);
  }
  vary /= *nsample;

  varyts = 0.0;
  avyts = 0.0;
  if (*testdat) {
    for (n = 0; n < ntest; ++n) {
      ntrue = n;
      varyts += ntrue * (yts[n] - avyts)*(yts[n] - avyts) / (ntrue + 1);
      avyts = (ntrue * avyts + yts[n]) / (ntrue + 1);
    }
    varyts /= ntest;
  }

  astr = 0.0;
  asd = 0.0;

  if (*iprox) {
    for (n = 0; n < *nsample; ++n) {
      for (k = 0; k < *nsample; ++k) {
	prox[k * *nsample + n] = 0.0;
	if (*oobprox) {
	  oobpair[k * *nsample + n] = 0;
	}
      }      
    }
    if (*testdat) {
      for (n = 0; n < ntest; ++n) {
	for (k = 0; k < ntest + *nsample; ++k) {
	  prox[k * ntest + n] = 0.0;
	}
      }
    }
  }

  for (m = 0; m < *mdim; ++m) {
    tgini[m] = 0.0;
  }
  
  if (*imp) {
    for (m = 0; m < *mdim; ++m) {
      errimp[m] = 0.0;
      errimp[m + *mdim] = 0.0;
      sse[m] = 0.0;
      /*
      for (n = 0; n < *nsample; ++n) {
	predimp[n + m * *nsample] = 0.0;
      }
      */
    }
  }

  if (*labelts) {
    for (n = 0; n < ntest; ++n) {
      ypred[n] = 0.0;
    }
  }

  GetRNGstate();

  /*************************************
   Start the loop over trees.
  *************************************/
  for (jb = 0; jb < *jbt; ++jb) {
    idx = (*keepf) ? jb * *nrnodes : 0;
    for (n = 0; n < *nsample; ++n) {
      jin[n] = 0;
      nodex[n] = 0;
    }

    for (n = 0; n < *nsample; ++n) {
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
    for (k = *nrnodes-1; k >= 0; --k) {
      if (nodestatus[k + idx] == 0) {
	ndbigtree[jb]--;
      }
      if (nodestatus[k + idx] == 2) {
	nodestatus[k + idx] = -1;
      }
    }

    for (n = 0; n < *nsample; ++n) {
      ytr[n] = 0.0;
    }
    
    F77_CALL(rtestreebag)(x, nsample, mdim, treemap + 2*idx, nodestatus + idx, 
			 nrnodes, ndbigtree + jb, ytr, upper + idx, 
			 avnode + idx, mbest + idx, cat, nodex);
    /* ytr is the prediction on OOB data by the current tree */
    /* yptr is the aggregated prediction by all trees grown so far */
    errb = 0.0;
    jout = 0;
    for (n = 0; n < *nsample; ++n) {
      if (jin[n] == 0) { 
	yptr[n] = (nout[n] * yptr[n] + ytr[n]) / (nout[n] + 1);
	nout[n] ++;
      }
      
      if (nout[n]) {
	jout ++;
      }
      errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
    }
    errb /= jout;

    ooberr = 0.0;
    for (n = 0; n < *nsample; ++n) {
      if (jin[n] == 0) {
	ooberr += (y[n]-avy) * ytr[n];
      }
    }

    /* Do simple linear regression of y on yhat for bias correction. */
    if (*biasCorr) {
      simpleLinReg(*nsample, yptr, y, coef, &errb, nout);
    }

    if (*testdat) {
      for (i = 0; i < ntest; ++i) {
	ytree[i] = 0.0;
	nodexts[i] = 0;
      }
      F77_CALL(rtestreebag)(xts, &ntest, mdim, treemap + 2*idx, 
			    nodestatus + idx, nrnodes, ndbigtree + jb,
			    ytree, upper + idx, avnode + idx, 
			    mbest + idx, cat, nodexts);
      /* ytree is the prediction for test data by the current tree */
      /* ypred is the aggregated prediction by all trees grown so far */
      errts = 0.0;
      for (n = 0; n < ntest; ++n) {
	ypred[n] = (jb * ypred[n] + ytree[n]) / (jb + 1);

      }

      if (*labelts) {
	for (n = 0; n < ntest; ++n) {
	  if (*biasCorr) {
	    resid = yts[n] - (coef[0] + coef[1]*ypred[n]);
	  } else {
	    resid = yts[n] - ypred[n];
	  }
	  errts += resid * resid;
	}
	errts /= ntest;
      }
    }

    if ((jb + 1) % *jprint == 0) {
      Rprintf("%d: ", jb + 1);
      if(*labelts == 1) Rprintf("MSE(Test)=%f  %%Var(y)=%7.2f  ", 
				errts, 100.0 * errts / varyts);
      Rprintf("MSE(OOB)=%f  %%Var(y)=%7.2f\n", errb, 100*errb/vary);
    }
    mse[jb] = errb;
    if(*labelts) {
      msets[jb] = errts;
    }

    /*  DO PROXIMITIES */
    if (*iprox) {
      for (n = 0; n < *nsample; ++n) {
	for (k = n + 1; k < *nsample; ++k) {
	  if (*oobprox) {
	    if (jin[n] == 0 && jin[k] == 0) {
	      oobpair[k * *nsample + n] ++;
	      oobpair[n * *nsample + k] ++;
	      if (nodex[k] == nodex[n]) {
		prox[k * *nsample + n] += 1.0;
		prox[n * *nsample + k] += 1.0;
	      }
	    }
	  } else {
	    if (nodex[k] == nodex[n]) {
	      prox[k * *nsample + n] += 1.0;
	      prox[n * *nsample + k] += 1.0;
	    }
	  }
	}
      }

      /* proximity for test data */
      if (*testdat) {
	for (n = 0; n < ntest; ++n) {
	  for (k = 0; k <= n; ++k) {
	    if (nodexts[k] == nodexts[n]) {
	      proxts[k * ntest + n] += 1.0;
	      proxts[n * ntest + k] = proxts[k * ntest + n];
	    }
	  }
	  for (k = 0; k < *nsample; ++k) {
	    if (nodexts[n] == nodex[k]) {
	      proxts[n + ntest * (k+ntest)] += 1.0; 
	    }
	  } 
	}
      } 
    }

    /* Variable importance */
    if (*imp) { 
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
	/* em = 0.0; */
	ooberrperm = 0.0;
	for (n = 0; n < *nsample; ++n) {
	  if (jin[n] == 0) {
	    ooberrperm += (y[n] - avy) * ytr[n];
	    /*
	    predimp[n + mr * *nsample] = (nout[n] * predimp[n + mr * *nsample]
					  + ytr[n]) / (nout[n] + 1);
	    em += (y[n] - predimp[n + mr* *nsample]) *
	      (y[n] - predimp[n + mr* *nsample]);
	    */
	  }
	}
	delta = (ooberr - ooberrperm) / *nsample;
	errimp[mr] += delta;
	sse[mr] += delta*delta;
	/* errimp[mr] = em / *nsample; */
      }
    }
  }
  PutRNGstate();

  /* end of tree iterations=======================================*/
  *rsq = 1.0 - errb/vary;

  if (*biasCorr) {
    for (n = 0; n < *nsample; ++n) {
      yptr[n] = coef[0] + coef[1]*yptr[n];
    }
    if (*testdat) {
      for (n = 0; n < ntest; ++n) {
	ypred[n] = coef[0] + coef[1]*ypred[n];
      }
    }
  }

  if (*iprox) {
    for (n = 0; n < *nsample; ++n) {
      for (k = n + 1; k < *nsample; ++k) {
	if (*oobprox) {
	  if (oobpair[k * *nsample + n] > 0) {
	    prox[k * *nsample + n] /= oobpair[k * *nsample + n];
	    prox[n * *nsample + k] /= oobpair[n * *nsample + k];
	  }
	} else {
	  prox[k * *nsample + n] /= oobpair[k * *nsample + n];
	  prox[n * *nsample + k] /= oobpair[n * *nsample + k];
	}
      }
      prox[n * *nsample + n] = 1.0;
    }
  }

  if (*imp) {
    for (m = 0; m < *mdim; ++m) {
      /* errimp[m] = 100 * ((errimp[m] / errb) - 1); */
      errimp[m] = errimp[m] / *jbt;
      sse[m] = sqrt( ((sse[m] / *jbt) - (errimp[m] * errimp[m])) / *jbt );
      errimp[m] = errimp[m] / sse[m];
      /*
	avimp(m)=avimp(m)/jbt
	sse(m)=(sse(m)/jbt)-avimp(m)*avimp(m)
	sse(m)=sqrt(sse(m)/jbt)
	zscore(m)=avimp(m)/sse(m)
      */
      errimp[m + *mdim] = tgini[m];
    }
  } else {
    for (m = 0; m < *mdim; ++m) {
      /* errimp[m] = tgini[m] / maximp; */
      errimp[m] = tgini[m];
    }
  }
  
}


void runrforest(double *xts, double *ypred, int *mdim, int *ntest, int *ntree, 
		int *ndbigtree, int *treemap, int *nodestatus, int *nrnodes, 
		double *upper, double *avnodes, int *mbest, int *cat,
		int *keepPred, double *allpred, int *iprox, double *proximity) 
{
  int i, j, k, n, idx, *nodex;
  double *ytree;

  ytree = (double *) S_alloc(*ntest, sizeof(double));
  nodex = (int *) S_alloc(*ntest, sizeof(int));

  if (*iprox) {
    for (i = 0; i < *ntest; ++i) {
      for (j = 0; j < *ntest; ++j) {
	proximity[i + j * *ntest] = 0.0;
	proximity[j + i * *ntest] = 0.0;
      }
    }
  }

  for(i = 0; i < *ntree; ++i) {

    if (*keepPred) {
      for (j = 0; j < *ntest; ++j) {
	allpred[j + i * *ntest] = 0.0;
      }
    }

    idx = i * *nrnodes;
    for(j = 0; j < *ntest; j++) ytree[j] = 0.0;
    F77_CALL(rtestreebag)(xts, ntest, mdim, treemap + 2*idx, nodestatus + idx,
			 nrnodes, ndbigtree + i, ytree, upper + idx, 
			 avnodes + idx, mbest + idx, cat, nodex);

    for(j = 0; j < *ntest; ++j) {
      ypred[j] += ytree[j];
    }

    if (*keepPred) {
      for(j = 0; j < *ntest; ++j) {
	allpred[j + i * *ntest] = ytree[j];
      }
    }

    /* if desired, do proximities for this round */
    if (*iprox) {
      for (n = 0; n < *ntest; ++n) {
	for (k = 0; k <= n; ++k) {
	  if(nodex[n] == nodex[k]) {
	    proximity[n + *ntest * k] += 1.0;
	    proximity[k + *ntest * n] = proximity[n + *ntest * k];
	  }
	}
      }
    }
  }
  for (i = 0; i < *ntest; ++i) {
    ypred[i] /= *ntree;
  }

  if (*iprox) {
    for (i = 0; i < *ntest; ++i) {
      for (j = 0; j < *ntest; ++j) {
	proximity[j + i * *ntest] /= *ntree;
      }
    }
  }
}

void simpleLinReg(int nsample, double *x, double *y, double *coef, 
		  double *mse, int *hasPred) {
  int i, nout = 0;
  double sxx=0.0, sxy=0.0, xbar=0.0, ybar=0.0;
  double dx = 0.0, dy = 0.0;

  for (i = 0; i < nsample; ++i) {
    if (hasPred[i]) {
      nout++;
      xbar += x[i];
      ybar += y[i];
    }
  }
  xbar /= nout;
  ybar /= nout;

  for (i = 0; i < nsample; ++i) {
    if (hasPred[i]) {
      dx = x[i] - xbar;
      dy = y[i] - ybar;
      sxx += dx * dx;
      sxy += dx * dy;
    }
  }
  coef[1] = sxy / sxx;
  coef[0] = ybar - coef[1] * xbar;
  
  *mse = 0.0;
  for (i = 0; i < nsample; ++i) {
    if (hasPred[i]) {
      dy = y[i] - (coef[0] + coef[1]*x[i]);
      *mse += dy * dy;
      /*      Rprintf("%f %f\n", y[i], dy); */
    }
  }
  *mse /= nout;

  return;
}
