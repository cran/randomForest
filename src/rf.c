/*****************************************************************
Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

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

#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include "rf.h"

void oob(int nsample, int nclass, int *jin, int *cl, int *jtr,int *jerr,
         int *counttr, int *out, double *errtr, int *jest, double *cutoff);

void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
		  int nclass, int nvote, double *errts,
		  int labelts, int *nclts, double *cutoff);

/*  Define the R RNG for use from Fortran. */
void F77_SUB(rrand)(double *r) { *r = unif_rand(); }

void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
	     int *sampsize, int *strata, int *Options, int *ntree, int *nvar,
	     int *ipi, double *classwt, double *cut, int *nodesize,
	     int *outcl, int *counttr, double *prox,
	     double *imprt, double *impsd, double *impmat, int *nrnodes,
	     int *ndbigtree, int *nodestatus, int *bestvar, int *treemap,
	     int *nodeclass, double *xbestsplit, double *errtr,
	     int *testdat, double *xts, int *clts, int *nts, double *countts,
	     int *outclts, int *labelts, double *proxts, double *errts,
             int *inbag) {
    /******************************************************************
     *  C wrapper for random forests:  get input from R and drive
     *  the Fortran routines.
     *
     *  Input:
     *
     *  x:        matrix of predictors (transposed!)
     *  dimx:     two integers: number of variables and number of cases
     *  cl:       class labels of the data
     *  ncl:      number of classes in the response
     *  cat:      integer vector of number of classes in the predictor;
     *            1=continuous
     * maxcat:    maximum of cat
     * Options:   7 integers: (0=no, 1=yes)
     *     add a second class (for unsupervised RF)?
     *         1: sampling from product of marginals
     *         2: sampling from product of uniforms
     *     assess variable importance?
     *     calculate proximity?
     *     calculate proximity based on OOB predictions?
     *     calculate outlying measure?
     *     how often to print output?
     *     keep the forest for future prediction?
     *  ntree:    number of trees
     *  nvar:     number of predictors to use for each split
     *  ipi:      0=use class proportion as prob.; 1=use supplied priors
     *  pi:       double vector of class priors
     *  nodesize: minimum node size: no node with fewer than ndsize
     *            cases will be split
     *
     *  Output:
     *
     *  outcl:    class predicted by RF
     *  counttr:  matrix of votes (transposed!)
     *  imprt:    matrix of variable importance measures
     *  impmat:   matrix of local variable importance measures
     *  prox:     matrix of proximity (if iprox=1)
     ******************************************************************/

    int nsample0, mdim, nclass, addClass, mtry, ntest, nsample, ndsize,
        mimp, nimp, near, nuse, noutall, nrightall, nrightimpall,
	keepInbag, nstrata;
    int jb, j, n, m, k, idxByNnode, idxByNsample, imp, localImp, iprox,
	oobprox, keepf, replace, stratify, trace, *nright,
	*nrightimp, *nout, *nclts, Ntree;

    int *out, *bestsplitnext, *bestsplit, *nodepop, *jin, *nodex,
	*nodexts, *nodestart, *ta, *ncase, *jerr, *varUsed,
	*jtr, *classFreq, *idmove, *jvr,
	*at, *a, *b, *mind, *nind, *jts, *oobpair;
    int **strata_idx, *strata_size, last, ktmp, nEmpty, ntry;

    double av=0.0, delta=0.0;

    double *tgini, *tx, *wl, *classpop, *tclasscat, *tclasspop, *win,
        *tp, *wr;

    addClass = Options[0];
    imp      = Options[1];
    localImp = Options[2];
    iprox    = Options[3];
    oobprox  = Options[4];
    trace    = Options[5];
    keepf    = Options[6];
    replace  = Options[7];
    stratify = Options[8];
    keepInbag = Options[9];
    mdim     = dimx[0];
    nsample0 = dimx[1];
    nclass   = (*ncl==1) ? 2 : *ncl;
    ndsize   = *nodesize;
    Ntree    = *ntree;
    mtry     = *nvar;
    ntest    = *nts;
    nsample = addClass ? (nsample0 + nsample0) : nsample0;
    mimp = imp ? mdim : 1;
    nimp = imp ? nsample : 1;
    near = iprox ? nsample0 : 1;
    if (trace == 0) trace = Ntree + 1;

    tgini =      (double *) S_alloc(mdim, sizeof(double));
    wl =         (double *) S_alloc(nclass, sizeof(double));
    wr =         (double *) S_alloc(nclass, sizeof(double));
    classpop =   (double *) S_alloc(nclass* *nrnodes, sizeof(double));
    tclasscat =  (double *) S_alloc(nclass*MAX_CAT, sizeof(double));
    tclasspop =  (double *) S_alloc(nclass, sizeof(double));
    tx =         (double *) S_alloc(nsample, sizeof(double));
    win =        (double *) S_alloc(nsample, sizeof(double));
    tp =         (double *) S_alloc(nsample, sizeof(double));

    out =           (int *) S_alloc(nsample, sizeof(int));
    bestsplitnext = (int *) S_alloc(*nrnodes, sizeof(int));
    bestsplit =     (int *) S_alloc(*nrnodes, sizeof(int));
    nodepop =       (int *) S_alloc(*nrnodes, sizeof(int));
    nodestart =     (int *) S_alloc(*nrnodes, sizeof(int));
    jin =           (int *) S_alloc(nsample, sizeof(int));
    nodex =         (int *) S_alloc(nsample, sizeof(int));
    nodexts =       (int *) S_alloc(ntest, sizeof(int));
    ta =            (int *) S_alloc(nsample, sizeof(int));
    ncase =         (int *) S_alloc(nsample, sizeof(int));
    jerr =          (int *) S_alloc(nsample, sizeof(int));
    varUsed =       (int *) S_alloc(mdim, sizeof(int));
    jtr =           (int *) S_alloc(nsample, sizeof(int));
    jvr =           (int *) S_alloc(nsample, sizeof(int));
    classFreq =     (int *) S_alloc(nclass, sizeof(int));
    jts =           (int *) S_alloc(ntest, sizeof(int));
    idmove =        (int *) S_alloc(nsample, sizeof(int));
    at =            (int *) S_alloc(mdim*nsample, sizeof(int));
    a =             (int *) S_alloc(mdim*nsample, sizeof(int));
    b =             (int *) S_alloc(mdim*nsample, sizeof(int));
    mind =          (int *) S_alloc(mdim, sizeof(int));
    nright =        (int *) S_alloc(nclass, sizeof(int));
    nrightimp =     (int *) S_alloc(nclass, sizeof(int));
    nout =          (int *) S_alloc(nclass, sizeof(int));
    if (oobprox) {
    	oobpair = (int *) S_alloc(near*near, sizeof(int));
    }

    /* Count number of cases in each class. */
    zeroInt(classFreq, nclass);
    for (n = 0; n < nsample; ++n) classFreq[cl[n] - 1] ++;
    /* Normalize class weights. */
    normClassWt(cl, nsample, nclass, *ipi, classwt, classFreq);

    if (stratify) {
	/* Count number of strata and frequency of each stratum. */
	nstrata = 0;
	for (n = 0; n < nsample0; ++n)
	    if (strata[n] > nstrata) nstrata = strata[n];
        /* Create the array of pointers, each pointing to a vector
	   of indices of where data of each stratum is. */
        strata_size = (int  *) S_alloc(nstrata, sizeof(int));
	for (n = 0; n < nsample0; ++n) {
	    strata_size[strata[n] - 1] ++;
	}
	strata_idx =  (int **) S_alloc(nstrata, sizeof(int *));
	for (n = 0; n < nstrata; ++n) {
	    strata_idx[n] = (int *) S_alloc(strata_size[n], sizeof(int));
	}
	zeroInt(strata_size, nstrata);
	for (n = 0; n < nsample0; ++n) {
	    strata_size[strata[n] - 1] ++;
	    strata_idx[strata[n] - 1][strata_size[strata[n] - 1] - 1] = n;
	}
    } else {
	nind = replace ? NULL : (int *) S_alloc(nsample, sizeof(int));
    }

    /*    INITIALIZE FOR RUN */
    if (*testdat) zeroDouble(countts, ntest * nclass);
    zeroInt(counttr, nclass * nsample);
    zeroInt(out, nsample);
    zeroDouble(tgini, mdim);
    zeroDouble(errtr, (nclass + 1) * Ntree);

    if (*labelts) {
	nclts  = (int *) S_alloc(nclass, sizeof(int));
	for (n = 0; n < ntest; ++n) nclts[clts[n]-1]++;
	zeroDouble(errts, (nclass + 1) * Ntree);
    }

    if (imp) {
        zeroDouble(imprt, (nclass+2) * mdim);
        zeroDouble(impsd, (nclass+1) * mdim);
	if (localImp) zeroDouble(impmat, nsample * mdim);
    }
    if (iprox) {
        zeroDouble(prox, nsample0 * nsample0);
        if (*testdat) zeroDouble(proxts, ntest * (ntest + nsample0));
    }
    makeA(x, mdim, nsample, cat, at, b);

    R_CheckUserInterrupt();

    /* Starting the main loop over number of trees. */
    GetRNGstate();
    if (trace <= Ntree) {
	/* Print header for running output. */
	Rprintf("ntree      OOB");
	for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	if (*labelts) {
	    Rprintf("|    Test");
	    for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	}
	Rprintf("\n");
    }
    idxByNnode = 0;
    idxByNsample = 0;
    for (jb = 0; jb < Ntree; jb++) {
        /* Do we need to simulate data for the second class? */
        if (addClass) createClass(x, nsample0, nsample, mdim);
		do {
			zeroInt(nodestatus + idxByNnode, *nrnodes);
			zeroInt(treemap + 2*idxByNnode, 2 * *nrnodes);
			zeroDouble(xbestsplit + idxByNnode, *nrnodes);
			zeroInt(nodeclass + idxByNnode, *nrnodes);
            zeroInt(varUsed, mdim);
            /* TODO: Put all sampling code into a function. */
            /* drawSample(sampsize, nsample, ); */
			if (stratify) {  /* stratified sampling */
				zeroInt(jin, nsample);
				zeroDouble(tclasspop, nclass);
				zeroDouble(win, nsample);
				if (replace) {  /* with replacement */
					for (n = 0; n < nstrata; ++n) {
						for (j = 0; j < sampsize[n]; ++j) {
							ktmp = (int) (unif_rand() * strata_size[n]);
							k = strata_idx[n][ktmp];
							tclasspop[cl[k] - 1] += classwt[cl[k] - 1];
							win[k] += classwt[cl[k] - 1];
							jin[k] = 1;
						}
					}
				} else { /* stratified sampling w/o replacement */
					/* re-initialize the index array */
					zeroInt(strata_size, nstrata);
					for (j = 0; j < nsample; ++j) {
						strata_size[strata[j] - 1] ++;
						strata_idx[strata[j] - 1][strata_size[strata[j] - 1] - 1] = j;
					}
					/* sampling without replacement */
					for (n = 0; n < nstrata; ++n) {
						last = strata_size[n] - 1;
						for (j = 0; j < sampsize[n]; ++j) {
							ktmp = (int) (unif_rand() * (last+1));
							k = strata_idx[n][ktmp];
                            swapInt(strata_idx[n][last], strata_idx[n][ktmp]);
							last--;
							tclasspop[cl[k] - 1] += classwt[cl[k]-1];
							win[k] += classwt[cl[k]-1];
							jin[k] = 1;
						}
					}
				}
			} else {  /* unstratified sampling */
				ntry = 0;
				do {
					nEmpty = 0;
					zeroInt(jin, nsample);
					zeroDouble(tclasspop, nclass);
					zeroDouble(win, nsample);
					if (replace) {
						for (n = 0; n < *sampsize; ++n) {
							k = unif_rand() * nsample;
							tclasspop[cl[k] - 1] += classwt[cl[k]-1];
							win[k] += classwt[cl[k]-1];
							jin[k] = 1;
						}
					} else {
						for (n = 0; n < nsample; ++n) nind[n] = n;
						last = nsample - 1;
						for (n = 0; n < *sampsize; ++n) {
							ktmp = (int) (unif_rand() * (last+1));
							k = nind[ktmp];
                            swapInt(nind[ktmp], nind[last]);
							last--;
							tclasspop[cl[k] - 1] += classwt[cl[k]-1];
							win[k] += classwt[cl[k]-1];
							jin[k] = 1;
						}
					}
					/* check if any class is missing in the sample */
					for (n = 0; n < nclass; ++n) {
						if (tclasspop[n] == 0.0) nEmpty++;
					}
					ntry++;
				} while (nclass - nEmpty < 2 && ntry <= 30);
				/* If there are still fewer than two classes in the data, throw an error. */
				if (nclass - nEmpty < 2) error("Still have fewer than two classes in the in-bag sample after 30 attempts.");
			}

            /* If need to keep indices of inbag data, do that here. */
            if (keepInbag) {
                for (n = 0; n < nsample0; ++n) {
                    inbag[n + idxByNsample] = jin[n];
                }
            }

			/* Copy the original a matrix back. */
			memcpy(a, at, sizeof(int) * mdim * nsample);
      	    modA(a, &nuse, nsample, mdim, cat, *maxcat, ncase, jin);

			F77_CALL(buildtree)(a, b, cl, cat, maxcat, &mdim, &nsample,
								&nclass,
								treemap + 2*idxByNnode, bestvar + idxByNnode,
								bestsplit, bestsplitnext, tgini,
								nodestatus + idxByNnode, nodepop,
								nodestart, classpop, tclasspop, tclasscat,
								ta, nrnodes, idmove, &ndsize, ncase,
								&mtry, varUsed, nodeclass + idxByNnode,
								ndbigtree + jb, win, wr, wl, &mdim,
								&nuse, mind);
			/* if the "tree" has only the root node, start over */
		} while (ndbigtree[jb] == 1);

		Xtranslate(x, mdim, *nrnodes, nsample, bestvar + idxByNnode,
				   bestsplit, bestsplitnext, xbestsplit + idxByNnode,
				   nodestatus + idxByNnode, cat, ndbigtree[jb]);

		/*  Get test set error */
		if (*testdat) {
            predictClassTree(xts, ntest, mdim, treemap + 2*idxByNnode,
                             nodestatus + idxByNnode, xbestsplit + idxByNnode,
                             bestvar + idxByNnode,
                             nodeclass + idxByNnode, ndbigtree[jb],
                             cat, nclass, jts, nodexts, *maxcat);
			TestSetError(countts, jts, clts, outclts, ntest, nclass, jb+1,
						 errts + jb*(nclass+1), *labelts, nclts, cut);
		}

		/*  Get out-of-bag predictions and errors. */
        predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
                         nodestatus + idxByNnode, xbestsplit + idxByNnode,
                         bestvar + idxByNnode,
                         nodeclass + idxByNnode, ndbigtree[jb],
                         cat, nclass, jtr, nodex, *maxcat);

		zeroInt(nout, nclass);
		noutall = 0;
		for (n = 0; n < nsample; ++n) {
			if (jin[n] == 0) {
				/* increment the OOB votes */
				counttr[n*nclass + jtr[n] - 1] ++;
				/* count number of times a case is OOB */
				out[n]++;
				/* count number of OOB cases in the current iteration.
				   nout[n] is the number of OOB cases for the n-th class.
				   noutall is the number of OOB cases overall. */
				nout[cl[n] - 1]++;
				noutall++;
			}
		}

        /* Compute out-of-bag error rate. */
		oob(nsample, nclass, jin, cl, jtr, jerr, counttr, out,
			errtr + jb*(nclass+1), outcl, cut);

		if ((jb+1) % trace == 0) {
			Rprintf("%5i: %6.2f%%", jb+1, 100.0*errtr[jb * (nclass+1)]);
			for (n = 1; n <= nclass; ++n) {
				Rprintf("%6.2f%%", 100.0 * errtr[n + jb * (nclass+1)]);
			}
			if (*labelts) {
				Rprintf("| ");
				for (n = 0; n <= nclass; ++n) {
					Rprintf("%6.2f%%", 100.0 * errts[n + jb * (nclass+1)]);
				}
			}
			Rprintf("\n");
#ifdef WIN32
			R_FlushConsole();
			R_ProcessEvents();
#endif
			R_CheckUserInterrupt();
		}

		/*  DO PROXIMITIES */
		if (iprox) {
            computeProximity(prox, oobprox, nodex, jin, oobpair, near);
			/* proximity for test data */
			if (*testdat) {
                computeProximity(proxts, 0, nodexts, jin, oobpair, ntest);
                /* Compute proximity between testset and training set. */
				for (n = 0; n < ntest; ++n) {
					for (k = 0; k < near; ++k) {
						if (nodexts[n] == nodex[k])
							proxts[n + ntest * (k+ntest)] += 1.0;
					}
				}
			}
		}

		/*  DO VARIABLE IMPORTANCE  */
		if (imp) {
			nrightall = 0;
			/* Count the number of correct prediction by the current tree
			   among the OOB samples, by class. */
			zeroInt(nright, nclass);
			for (n = 0; n < nsample; ++n) {
       	        /* out-of-bag and predicted correctly: */
				if (jin[n] == 0 && jtr[n] == cl[n]) {
					nright[cl[n] - 1]++;
					nrightall++;
				}
			}
			for (m = 0; m < mdim; ++m) {
				if (varUsed[m]) {
					nrightimpall = 0;
					zeroInt(nrightimp, nclass);
					for (n = 0; n < nsample; ++n) tx[n] = x[m + n*mdim];
					/* Permute the m-th variable. */
                    permuteOOB(m, x, jin, nsample, mdim);
					/* Predict the modified data using the current tree. */
                    predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
                                     nodestatus + idxByNnode,
                                     xbestsplit + idxByNnode,
                                     bestvar + idxByNnode,
                                     nodeclass + idxByNnode, ndbigtree[jb],
                                     cat, nclass, jvr, nodex, *maxcat);
					/* Count how often correct predictions are made with
					   the modified data. */
					for (n = 0; n < nsample; n++) {
						/* Restore the original data for that variable. */
						x[m + n*mdim] = tx[n];
						if (jin[n] == 0) {
							if (jvr[n] == cl[n]) {
								nrightimp[cl[n] - 1]++;
								nrightimpall++;
							}
							if (localImp && jvr[n] != jtr[n]) {
								if (cl[n] == jvr[n]) {
									impmat[m + n*mdim] -= 1.0;
								} else {
									impmat[m + n*mdim] += 1.0;
								}
							}
						}
					}
					/* Accumulate decrease in proportions of correct
					   predictions. */
					/* class-specific measures first: */
					for (n = 0; n < nclass; ++n) {
						if (nout[n] > 0) {
							delta = ((double) (nright[n] - nrightimp[n])) / nout[n];
							imprt[m + n*mdim] += delta;
							impsd[m + n*mdim] += delta * delta;
						}
					}
					/* overall measure, across all classes: */
					if (noutall > 0) {
						delta = ((double)(nrightall - nrightimpall)) / noutall;
						imprt[m + nclass*mdim] += delta;
						impsd[m + nclass*mdim] += delta * delta;
					}
				}
			}
		}

		R_CheckUserInterrupt();
#ifdef WIN32
		R_ProcessEvents();
#endif
        if (keepf) idxByNnode += *nrnodes;
        if (keepInbag) idxByNsample += nsample0;
    }
    PutRNGstate();

    /*  Final processing of variable importance. */
    for (m = 0; m < mdim; m++) tgini[m] /= Ntree;
    if (imp) {
		for (m = 0; m < mdim; ++m) {
			if (localImp) { /* casewise measures */
				for (n = 0; n < nsample; ++n) impmat[m + n*mdim] /= out[n];
			}
			/* class-specific measures */
			for (k = 0; k < nclass; ++k) {
				av = imprt[m + k*mdim] / Ntree;
				impsd[m + k*mdim] =
                    sqrt(((impsd[m + k*mdim] / Ntree) - av*av) / Ntree);
				imprt[m + k*mdim] = av;
				/* imprt[m + k*mdim] = (se <= 0.0) ? -1000.0 - av : av / se; */
			}
			/* overall measures */
			av = imprt[m + nclass*mdim] / Ntree;
			impsd[m + nclass*mdim] =
                sqrt(((impsd[m + nclass*mdim] / Ntree) - av*av) / Ntree);
			imprt[m + nclass*mdim] = av;
			imprt[m + (nclass+1)*mdim] = tgini[m];
		}
    } else {
		for (m = 0; m < mdim; ++m) imprt[m] = tgini[m];
    }

    /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/
    if (iprox) {
		for (n = 0; n < near; ++n) {
			for (k = n + 1; k < near; ++k) {
                prox[near*k + n] /= oobprox ?
                    (oobpair[near*k + n] > 0 ? oobpair[near*k + n] : 1) :
                    Ntree;
				prox[near*n + k] = prox[near*k + n];
			}
			prox[near*n + n] = 1.0;
		}
		if (*testdat) {
			for (n = 0; n < ntest; ++n) {
				for (k = 0; k < ntest + nsample; ++k)
					proxts[ntest*k + n] /= Ntree;
				proxts[ntest * n + n] = 1.0;
			}
		}
    }
}


void classForest(int *mdim, int *ntest, int *nclass, int *maxcat,
                 int *nrnodes, int *ntree, double *x, double *xbestsplit,
                 double *pid, double *cutoff, double *countts, int *treemap,
                 int *nodestatus, int *cat, int *nodeclass, int *jts,
                 int *jet, int *bestvar, int *node, int *treeSize,
                 int *keepPred, int *prox, double *proxMat, int *nodes) {
    int j, n, n1, n2, idxNodes, offset1, offset2, *junk, ntie;
    double crit, cmax;

    zeroDouble(countts, *nclass * *ntest);
    idxNodes = 0;
    offset1 = 0;
    offset2 = 0;
    junk = NULL;

    for (j = 0; j < *ntree; ++j) {
		/* predict by the j-th tree */
        predictClassTree(x, *ntest, *mdim, treemap + 2*idxNodes,
			 nodestatus + idxNodes, xbestsplit + idxNodes,
			 bestvar + idxNodes, nodeclass + idxNodes,
			 treeSize[j], cat, *nclass,
			 jts + offset1, node + offset2, *maxcat);
		/* accumulate votes: */
		for (n = 0; n < *ntest; ++n) {
			countts[jts[n + offset1] - 1 + n * *nclass] += 1.0;
		}

		/* if desired, do proximities for this round */
		if (*prox) computeProximity(proxMat, 0, node + offset2, junk, junk,
									*ntest);
		idxNodes += *nrnodes;
		if (*keepPred) offset1 += *ntest;
		if (*nodes)    offset2 += *ntest;
    }

    /* Aggregated prediction is the class with the maximum votes/cutoff */
    for (n = 0; n < *ntest; ++n) {
		cmax = 0.0;
		ntie = 1;
		for (j = 0; j < *nclass; ++j) {
			crit = (countts[j + n * *nclass] / *ntree) / cutoff[j];
			if (crit > cmax) {
				jet[n] = j + 1;
				cmax = crit;
				ntie = 1;
			}
			/* Break ties at random: */
			if (crit == cmax) {
				if (unif_rand() < 1.0 / ntie) jet[n] = j + 1;
				ntie++;
			}
		}
    }

    /* if proximities requested, do the final adjustment
       (division by number of trees) */
    if (*prox) {
		for (n1 = 0; n1 < *ntest; ++n1) {
			for (n2 = n1 + 1; n2 < *ntest; ++n2) {
				proxMat[n1 + n2 * *ntest] /= *ntree;
				proxMat[n2 + n1 * *ntest] = proxMat[n1 + n2 * *ntest];
			}
			proxMat[n1 + n1 * *ntest] = 1.0;
		}
    }
}

/*
  Modified by A. Liaw 1/10/2003 (Deal with cutoff)
  Re-written in C by A. Liaw 3/08/2004
*/
void oob(int nsample, int nclass, int *jin, int *cl, int *jtr,int *jerr,
	 int *counttr, int *out, double *errtr, int *jest,
	 double *cutoff) {
    int j, n, noob, *noobcl, ntie;
    double qq, smax, smaxtr;

    noobcl  = (int *) S_alloc(nclass, sizeof(int));
    zeroInt(jerr, nsample);
    zeroDouble(errtr, nclass+1);

    noob = 0;
    for (n = 0; n < nsample; ++n) {
	if (out[n]) {
	    noob++;
	    noobcl[cl[n]-1]++;
	    smax = 0.0;
	    smaxtr = 0.0;
	    ntie = 1;
	    for (j = 0; j < nclass; ++j) {
	    	qq = (((double) counttr[j + n*nclass]) / out[n]) / cutoff[j];
	    	if (j+1 != cl[n]) smax = (qq > smax) ? qq : smax;
	    	/* if vote / cutoff is larger than current max, re-set max and
		   	   change predicted class to the current class */
	    	if (qq > smaxtr) {
	    		smaxtr = qq;
	    		jest[n] = j+1;
	    		ntie = 1;
	    	}
	    	/* break tie at random */
	    	if (qq == smaxtr) {
	    		if (unif_rand() < 1.0 / ntie) {
	    			smaxtr = qq;
	    			jest[n] = j+1;
	    		}
	    		ntie++;
	    	}
	    }
	    if (jest[n] != cl[n]) {
		errtr[cl[n]] += 1.0;
		errtr[0] += 1.0;
		jerr[n] = 1;
	    }
	}
    }
    errtr[0] /= noob;
    for (n = 1; n <= nclass; ++n) errtr[n] /= noobcl[n-1];
}


void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
		  int nclass, int nvote, double *errts,
		  int labelts, int *nclts, double *cutoff) {
    int j, n, ntie;
    double cmax, crit;

    for (n = 0; n < ntest; ++n) countts[jts[n]-1 + n*nclass] += 1.0;

    /*  Prediction is the class with the maximum votes */
    for (n = 0; n < ntest; ++n) {
		cmax=0.0;
		ntie = 1;
		for (j = 0; j < nclass; ++j) {
			crit = (countts[j + n*nclass] / nvote) / cutoff[j];
			if (crit > cmax) {
				jet[n] = j+1;
				cmax = crit;
				ntie = 1;
			}
			/*  Break ties at random: */
			if (crit == cmax) {
				if (unif_rand() < 1.0 / ntie) {
					jet[n] = j+1;
					cmax = crit;
				}
				ntie++;
			}
		}
    }
    if (labelts) {
        zeroDouble(errts, nclass + 1);
		for (n = 0; n < ntest; ++n) {
			if (jet[n] != clts[n]) {
				errts[0] += 1.0;
				errts[clts[n]] += 1.0;
			}
		}
		errts[0] /= ntest;
		for (n = 1; n <= nclass; ++n) errts[n] /= nclts[n-1];
    }
}
