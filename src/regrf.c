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

#include <R.h>
#include "rf.h"

void simpleLinReg(int nsample, double *x, double *y, double *coef, 
		  double *mse, int *hasPred);


void regRF(double *x, double *y, int *nsample, int *mdim, int *sampsize,
	   int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp, 
	   int *cat, int *jprint, int *doProx, int *oobprox, int *biasCorr, 
	   double *yptr, double *errimp, double *impmat, double *impSD,
	   double *prox, int *ndbigtree, int *nodestatus, int *lDaughter, 
           int *rDaughter, double *avnode, int *mbest, double *upper, 
           double *mse, int *keepf, int *replace, int *testdat, double *xts, 
           int *nts, double *yts, int *labelts, double *yTestPred, 
           double *proxts, double *msets, double *coef, int *nout) {
    /*************************************************************************
   Input:
   mdim=number of variables in data set
   nsample=number of cases

   nthsize=number of cases in a node below which the tree will not split,
   setting nthsize=5 generally gives good results.

   nTree=number of trees in run.  200-500 gives pretty good results

   mtry=number of variables to pick to split on at each node.  mdim/3
   seems to give genrally good performance, but it can be 
   altered up or down

   imp=1 turns on variable importance.  This is computed for the
   mth variable as the percent rise in the test set mean sum-of-
   squared errors when the mth variable is randomly permuted.

  *************************************************************************/
  
    double errts = 0.0, averrb, meanY, meanYts, varY, varYts, xrand, 
	errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;
    
    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini; 
    
    int k, m, mr, n, jout, j, idx, ntest, last, ktmp;
    int *oobpair, varImp, localImp, *ind, *indts, *varUsed;
    
    int *in, *nind, *nodex, *nodexts;
    
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    
    if (*jprint == 0) *jprint = *nTree + 1;
    
    yb         = (double *) S_alloc(*sampsize, sizeof(double));
    xb         = (double *) S_alloc(*mdim * *sampsize, sizeof(double));
    ytr        = (double *) S_alloc(*nsample, sizeof(double));
    xtmp       = (double *) S_alloc(*nsample, sizeof(double));
    resOOB     = (double *) S_alloc(*nsample, sizeof(double));
    
    in        = (int *) S_alloc(*nsample, sizeof(int));
    ind        = (int *) S_alloc(*nsample, sizeof(int));
    nodex      = (int *) S_alloc(*nsample, sizeof(int));
    varUsed    = (int *) S_alloc(*mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(*nsample, sizeof(int));
    
    if (*testdat) {
	ytree      = (double *) S_alloc(ntest, sizeof(double));
	nodexts    = (int *) S_alloc(ntest, sizeof(int));
	indts      = (int *) S_alloc(ntest, sizeof(int));
    }
    oobpair = (*doProx && *oobprox) ? 
	(int *) S_alloc(*nsample * *nsample, sizeof(int)) : NULL;
    
    /* If variable importance is requested, tgini points to the second
       "column" of errimp, otherwise it's just the same as errimp. */
    tgini = varImp ? errimp + *mdim : errimp;

    averrb = 0.0;
    meanY = 0.0;
    varY = 0.0;
    
    zeroDouble(yptr, *nsample);
    zeroInt(nout, *nsample);
    zeroInt(ind, *nsample);
    for (n = 0; n < *nsample; ++n) {
	varY += n * (y[n] - meanY)*(y[n] - meanY) / (n + 1);
	meanY = (n * meanY + y[n]) / (n + 1);
    }
    varY /= *nsample;
    
    varYts = 0.0;
    meanYts = 0.0;
    if (*testdat) {
	for (n = 0; n < ntest; ++n) {
	    varYts += n * (yts[n] - meanYts)*(yts[n] - meanYts) / (n + 1);
	    meanYts = (n * meanYts + yts[n]) / (n + 1);
	}
	varYts /= ntest;
    }
    
    if (*doProx) {
	if (*testdat) {
	    zeroDouble(prox, *nsample * (*nsample + ntest));
	} else {
	    zeroDouble(prox, *nsample * *nsample);
	}
    }
    
    if (varImp) {
        zeroDouble(errimp, *mdim * 2);
	if (localImp) zeroDouble(impmat, *nsample * *mdim);
    } else {
        zeroDouble(errimp, *mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);
    
    /* print header for running output */
    if (*jprint <= *nTree) {
	Rprintf("     |      Out-of-bag   ");
	if (*testdat) {
	    Rprintf("|       Test set    ");
	}
	Rprintf("|\n");
	Rprintf("Tree |      MSE  %%Var(y) ");
	if (*testdat) {
	    Rprintf("|      MSE  %%Var(y) ");
	}
	Rprintf("|\n");
    }
    GetRNGstate();
    /*************************************
     * Start the loop over trees.
     *************************************/
    for (j = 0; j < *nTree; ++j) {
	idx = (*keepf) ? j * *nrnodes : 0;
	zeroInt(in, *nsample);
	zeroInt(nodex, *nsample);
        zeroInt(varUsed, *mdim);
	if (*replace) {
	    for (n = 0; n < *sampsize; ++n) {
		xrand = unif_rand();
		k = xrand * *nsample;
		in[k] = 1;
		yb[n] = y[k];
		for(m = 0; m < *mdim; ++m) {
		    xb[m + n * *mdim] = x[m + k * *mdim];
		}
	    }
	} else {
	    for (n = 0; n < *nsample; ++n) nind[n] = n;
	    last = *nsample - 1;
	    for (n = 0; n < *sampsize; ++n) {
		ktmp = (int) (unif_rand() * (last+1));
		k = nind[ktmp];
		nind[ktmp] = nind[last];
		nind[last] = k;
		last--;
		in[k] = 1;
		yb[n] = y[k];
		for(m = 0; m < *mdim; ++m) {
		    xb[m + n * *mdim] = x[m + k * *mdim];
		}
	    }
	}
	
	regTree(xb, yb, *mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodestatus + idx, *nrnodes, 
                *nthsize, *mtry, mbest + idx, cat, tgini, varUsed);
	ndbigtree[j] = *nrnodes;
	for (k = *nrnodes-1; k >= 0; --k) {
	    if (nodestatus[k + idx] == 0) {
		ndbigtree[j]--;
	    }
	    if (nodestatus[k + idx] == NODE_TOSPLIT) {
		nodestatus[k + idx] = NODE_TERMINAL;
	    }
	}

	zeroDouble(ytr, *nsample);	
	predictRegTree(x, *nsample, *mdim, ind, lDaughter + idx, 
                       rDaughter + idx, nodestatus + idx, ytr, upper + idx, 
                       avnode + idx, mbest + idx, cat, nodex);

	/* ytr is the prediction on OOB data by the current tree */
	/* yptr is the aggregated prediction by all trees grown so far */
	errb = 0.0;
	jout = 0;
	for (n = 0; n < *nsample; ++n) {
	    if (in[n] == 0) {
		nout[n]++;
		yptr[n] = ((nout[n]-1) * yptr[n] + ytr[n]) / nout[n];
		resOOB[n] = ytr[n] - y[n];
	    }
            if (nout[n]) jout ++;
	    errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
	}
	errb /= jout;
	
	ooberr = 0.0;
	for (n = 0; n < *nsample; ++n) {
	    if (in[n] == 0) ooberr += (y[n]-meanY) * ytr[n];
	}
	
	/* Do simple linear regression of y on yhat for bias correction. */
	if (*biasCorr) simpleLinReg(*nsample, yptr, y, coef, &errb, nout);
	
	if (*testdat) {
	    zeroDouble(ytree, ntest);
	    zeroInt(nodexts, ntest);
	    predictRegTree(xts, ntest, *mdim, indts, lDaughter + idx,
			   rDaughter + idx, nodestatus + idx, ytree, 
                           upper + idx, avnode + idx, 
			   mbest + idx, cat, nodexts);
	    /* ytree is the prediction for test data by the current tree */
	    /* yTestPred is the aggregated prediction by all trees grown so far */
	    errts = 0.0;
	    for (n = 0; n < ntest; ++n) {
		yTestPred[n] = (j * yTestPred[n] + ytree[n]) / (j + 1);
	    }
	    
	    if (*labelts) {
		for (n = 0; n < ntest; ++n) {
		    resid = *biasCorr ? yts[n] - (coef[0] + coef[1]*yTestPred[n]) :
                        yts[n] - yTestPred[n];
		    errts += resid * resid;
		}
		errts /= ntest;
	    }
	}
	
	/* Print running output. */
	if ((j + 1) % *jprint == 0) {
	    Rprintf("%4d |", j + 1);
	    Rprintf(" %8.4g %8.2f ", errb, 100 * errb / varY);
	    if(*labelts == 1) Rprintf("| %8.4g %8.2f ", 
				      errts, 100.0 * errts / varYts);
	    Rprintf("|\n");
	}
	mse[j] = errb;
	if(*labelts) {
	    msets[j] = errts;
	}
	
	/*  DO PROXIMITIES */
	if (*doProx) {
	    for (n = 0; n < *nsample; ++n) {
		for (k = n + 1; k < *nsample; ++k) {
		    if (*oobprox && in[n] == 0 && in[k] == 0) {
                        oobpair[k * *nsample + n] ++;
                        oobpair[n * *nsample + k] ++;
                    }
                    if (nodex[k] == nodex[n]) {
                        prox[k * *nsample + n] += 1.0;
                        prox[n * *nsample + k] += 1.0;
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
	if (varImp) { 
	    for (mr = 0; mr < *mdim; ++mr) {
                if (varUsed[mr]) { /* Go ahead if the variable is used */
                    /* make a copy of the m-th variable into xtmp */
                    for (n = 0; n < *nsample; ++n) xtmp[n] = x[mr + n * *mdim];
                    permuteOOB(mr, x, in, *nsample, *mdim);
                    predictRegTree(x, *nsample, *mdim, in, lDaughter + idx, 
                                   rDaughter + idx, nodestatus + idx, ytr, 
                                   upper + idx, avnode + idx, mbest + idx, 
                                   cat, nodex);
                    ooberrperm = 0.0;
                    for (n = 0; n < *nsample; ++n) {
                        x[mr + n * *mdim] = xtmp[n]; /* copy original data back */
                        if (in[n] == 0) {
                            ooberrperm += (y[n] - meanY) * ytr[n];
                            if (localImp) {
                                impmat[mr + n * *mdim] += 
                                    (ytr[n] - y[n]) * (ytr[n] - y[n]) - 
                                    resOOB[n] * resOOB[n];
                            }
                        }
                    }
                    delta = (ooberr - ooberrperm) / *nsample;
                    errimp[mr] += delta;
                    impSD[mr] += delta * delta;
                }
            }
        }
    }
    PutRNGstate();
    /* end of tree iterations=======================================*/
    
    if (*doProx) {
	for (n = 0; n < *nsample; ++n) {
	    for (k = n + 1; k < *nsample; ++k) {
		if (*oobprox) {
		    if (oobpair[k * *nsample + n]) {
			prox[k * *nsample + n] /= oobpair[k * *nsample + n];
			prox[n * *nsample + k] /= oobpair[n * *nsample + k];
		    }
		} else {
		    prox[k * *nsample + n] /= *nTree;
		    prox[n * *nsample + k] /= *nTree;
		}
	    }
	    prox[n * *nsample + n] = 1.0;
	}
        if (*testdat) {
            for (n = 0; n < ntest; ++n) {
                for (k = 0; k <= n; ++k) {
                    if (nodexts[k] == nodexts[n]) {
                        proxts[k * ntest + n] /= *nTree;
                        proxts[n * ntest + k] /= *nTree;
                    }
                }
                for (k = 0; k < *nsample; ++k) {
                    if (nodexts[n] == nodex[k]) {
                        proxts[n + ntest * (k+ntest)] /= *nTree; 
                    }
                } 
            }
        }
    }
    
    if (varImp) {
	for (m = 0; m < *mdim; ++m) {
	    errimp[m] = errimp[m] / *nTree;
	    impSD[m] = sqrt( ((impSD[m] / *nTree) - 
			      (errimp[m] * errimp[m])) / *nTree );
	    if (localImp) {
                for (n = 0; n < *nsample; ++n) {
                    impmat[m + n * *mdim] /= nout[n];
                }
	    }
        }
    }
    for (m = 0; m < *mdim; ++m) tgini[m] /= *nTree;
}

/*----------------------------------------------------------------------*/
void runrforest(double *xts, double *ypred, int *mdim, int *ntest, 
		int *ntree, int *ndbigtree, int *lDaughter, int *rDaughter, 
		int *nodestatus, int *nrnodes, double *upper, 
		double *avnodes, int *mbest, int *cat, int *keepPred, 
		double *allpred, int *doProx, double *proximity) {
    int i, j, k, n, idx, *nodex, *ind;
    double *ytree;

    ytree = (double *) S_alloc(*ntest, sizeof(double));
    nodex = (int *) S_alloc(*ntest, sizeof(int));
    ind   = (int *) S_alloc(*ntest, sizeof(int));

    if (*doProx) zeroDouble(proximity, *ntest * *ntest);
    if (*keepPred) zeroDouble(allpred, *ntest * *ntree);

    for(i = 0; i < *ntree; ++i) {
	idx = i * *nrnodes;
	zeroDouble(ytree, *ntest);
	predictRegTree(xts, *ntest, *mdim, ind, lDaughter + idx, 
                       rDaughter + idx,
		       nodestatus + idx, ytree, upper + idx, avnodes + idx, 
                       mbest + idx, cat, nodex);
	for(j = 0; j < *ntest; ++j) {
	    ypred[j] += ytree[j];
	}
	if (*keepPred) {
	    for(j = 0; j < *ntest; ++j) {
		allpred[j + i * *ntest] = ytree[j];
	    }
	}

	/* if desired, do proximities for this round */
	if (*doProx) {
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
    for (i = 0; i < *ntest; ++i) ypred[i] /= *ntree;
    if (*doProx) {
	for (i = 0; i < *ntest; ++i) {
	    for (j = 0; j < *ntest; ++j) {
		proximity[j + i * *ntest] /= *ntree;
	    }
	}
    }
}

void simpleLinReg(int nsample, double *x, double *y, double *coef, 
		  double *mse, int *hasPred) {
/* Compute simple linear regression of y on x, returning the coefficients,
   the average squared residual, and the predicted values (overwriting y). */
    int i, nout = 0;
    double sxx=0.0, sxy=0.0, xbar=0.0, ybar=0.0;
    double dx = 0.0, dy = 0.0, py=0.0;

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
            py = coef[0] + coef[1] * x[i];
	    dy = y[i] - py;
	    *mse += dy * dy;
            y[i] = py;
	}
    }
    *mse /= nout;
    return;
}
