/*****************************************************************
Copyright (C) 2001-4 Leo Breiman, Adele Cutler and Merck & Co., Inc.
  
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
#include "rf.h"

void oob(int nsample, int nclass, int *jin, int *cl, int *jtr,int *jerr,
         int *counttr, int *out, double *errtr, int *jest, double *cutoff);

void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
		  int nclass, int nvote, double *errts, double *pid, 
		  int labelts, int *nclts, double *cutoff);

/*  Define the R RNG for use from Fortran. */
void F77_SUB(rrand)(double *r) { *r = unif_rand(); }

void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat, 
	     int *sampsize, int *Options, int *ntree, int *nvar,
	     int *ipi, double *pi, double *cut, int *nodesize, 
	     int *outcl, int *counttr, double *prox, 
	     double *imprt, double *impsd, double *impmat, int *nrnodes, 
	     int *ndbigtree, int *nodestatus, int *bestvar, int *treemap, 
	     int *nodeclass, double *xbestsplit, double *pid, double *errtr, 
	     int *testdat, double *xts, int *clts, int *nts, double *countts,
	     int *outclts, int *labelts, double *proxts, double *errts) {
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

    int nsample0, mdim, nclass, addClass, jbt, mtry, ntest, nsample, ndsize,
        mimp, nimp, near, nuse, noutall, nrightall, nrightimpall;
    int jb, j, n, m, k, arrayindex, imp, localImp, iprox, 
	oobprox, keepf, replace, stratify, trace, *nright, 
	*nrightimp, *nout, *nclts;

    int *out, *bestsplitnext, *bestsplit, *zeroes,
	*nodepop, *parent, *jin, *nodex,
	*nodexts, *nodestart, *ta, *ncase, *jerr, *varUsed,
	*jtr, *nc, *idmove, *jvr,
	*at, *a, *b, *cbestsplit, *mind, *nind, *jts, *oobpair;
    int **class_idx, *class_size, last, tmp, ktmp, anyEmpty, ntry;

    double av=0.0;
    
    double *tgini, *tx, *wl, *classpop, *tclasscat, *tclasspop, *win,
        *tp, *wc, *wr, *wtt, *iw;

    addClass = Options[0];
    imp      = Options[1];
    localImp = Options[2];
    iprox    = Options[3];
    oobprox  = Options[4];
    trace    = Options[5]; 
    keepf    = Options[6];
    replace  = Options[7];
    stratify = Options[8];
    mdim     = dimx[0];
    nsample0 = dimx[1];
    nclass   = (*ncl==1) ? 2 : *ncl;
    ndsize   = *nodesize;
    jbt      = *ntree;
    mtry     = *nvar;
    ntest    = *nts; 
    nsample = addClass ? (nsample0 + nsample0) : nsample0;
    mimp = imp ? mdim : 1;
    nimp = imp ? nsample : 1;
    near = iprox ? nsample0 : 1;
    if (trace == 0) trace = jbt + 1;


    tgini =      (double *) S_alloc(mdim, sizeof(double));
    wl =         (double *) S_alloc(nclass, sizeof(double));
    wc =         (double *) S_alloc(nclass, sizeof(double));
    wr =         (double *) S_alloc(nclass, sizeof(double));
    classpop =   (double *) S_alloc(nclass* *nrnodes, sizeof(double));
    tclasscat =  (double *) S_alloc(nclass*32, sizeof(double));
    tclasspop =  (double *) S_alloc(nclass, sizeof(double));
    tx =         (double *) S_alloc(nsample, sizeof(double));
    win =        (double *) S_alloc(nsample, sizeof(double));
    tp =         (double *) S_alloc(nsample, sizeof(double));
    wtt =        (double *) S_alloc(nsample, sizeof(double));
    iw =         (double *) S_alloc(nsample, sizeof(double));

    out =           (int *) S_alloc(nsample, sizeof(int));
    bestsplitnext = (int *) S_alloc(*nrnodes, sizeof(int));
    bestsplit =     (int *) S_alloc(*nrnodes, sizeof(int));
    nodepop =       (int *) S_alloc(*nrnodes, sizeof(int));
    parent =        (int *) S_alloc(*nrnodes, sizeof(int));
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
    nc =            (int *) S_alloc(nclass, sizeof(int));
    jts =           (int *) S_alloc(ntest, sizeof(int));
    zeroes =        (int *) S_alloc(ntest, sizeof(int));
    idmove =        (int *) S_alloc(nsample, sizeof(int));
    at =            (int *) S_alloc(mdim*nsample, sizeof(int));
    a =             (int *) S_alloc(mdim*nsample, sizeof(int));
    b =             (int *) S_alloc(mdim*nsample, sizeof(int));
    cbestsplit =    (int *) S_alloc(*maxcat * *nrnodes, sizeof(int));
    mind =          (int *) S_alloc(mdim, sizeof(int));
    nright =        (int *) S_alloc(nclass, sizeof(int));
    nrightimp =     (int *) S_alloc(nclass, sizeof(int));
    nout =          (int *) S_alloc(nclass, sizeof(int));
    if (oobprox) {
	oobpair = (int *) S_alloc(near*near, sizeof(int));
    }

    prepare(cl, nsample, nclass, *ipi, pi, pid, nc, wtt);    

    if (stratify) {
        /* Create the array of pointers, each pointing to a vector 
	   of indices of where data of each class is. */
        class_size = (int  *) S_alloc(nclass, sizeof(int));
	class_idx =  (int **) S_alloc(nclass, sizeof(int *));
	for (n = 0; n < nclass; ++n) {
	    class_idx[n] = (int *) S_alloc(nc[n], sizeof(int));
	}
	for (n = 0; n < nsample; ++n) {
	    class_size[cl[n] - 1] ++;
	    class_idx[cl[n]-1][class_size[cl[n]-1] - 1] = n;
	}
    } else {
	if (replace) {
	    nind = NULL;
	} else {
	    nind = (int *) S_alloc(nsample, sizeof(int));
	}
    }

    /*    INITIALIZE FOR RUN */
    if (*testdat) {
        zeroDouble(countts, ntest * nclass);
        zeroInt(zeroes, ntest);
    }
    zeroInt(counttr, nclass * nsample);
    zeroInt(out, nsample);
    zeroDouble(tgini, mdim);
    zeroDouble(errtr, (nclass + 1) * jbt);

    if (*labelts) {
	nclts  = (int *) S_alloc(nclass, sizeof(int));
	for (n = 0; n < ntest; ++n) nclts[clts[n]-1]++;
	zeroDouble(errts, (nclass + 1) * jbt);
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

    /*   START RUN   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    GetRNGstate(); 
    if (trace <= jbt) {
	/* Print header for running output. */
	Rprintf("ntree      OOB");
	for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	if (*labelts) {
	    Rprintf("|    Test");
	    for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	}
	Rprintf("\n");
    }

    for(jb = 0; jb < jbt; jb++) {
        /* Do we need to simulate data for the second class? */
        /* SET UP DATA TO ADD A CLASS++++++++++++++++++++++++++++++ */
        if (addClass) createClass(x, nsample0, nsample, mdim);

	arrayindex = keepf ? jb * *nrnodes : 0;
	do {
	    zeroInt(nodestatus + arrayindex, *nrnodes);
	    zeroInt(treemap + 2*arrayindex, 2 * *nrnodes);
	    zeroDouble(xbestsplit + arrayindex, *nrnodes);
	    zeroInt(nodeclass + arrayindex, *nrnodes);
	    zeroInt(jin, nsample);
            zeroInt(varUsed, mdim);
	    zeroDouble(tclasspop, nclass);
	    zeroDouble(win, nsample);
      
	    if (stratify) {
		/* stratified sampling by class */
		if (replace) {
		    for (n = 0; n < nclass; ++n) {
			for (j = 0; j < sampsize[n]; ++j) {
			    ktmp = (int) (unif_rand() * nc[n]);
			    k = class_idx[n][ktmp];
			    tclasspop[cl[k] - 1] += wtt[k];
			    win[k] += wtt[k];
			    jin[k] = 1;
			}
		    }
		} else { /* stratified sampling w/o replacement */
		    /* re-initialize the index array */
  		    zeroInt(class_size, nclass);
		    for (j = 0; j < nsample; ++j) {
			class_size[cl[j] - 1] ++;
			class_idx[cl[j]-1][class_size[cl[j]-1] - 1] = j;
		    }
		    /* sampling without replacement */
		    for (n = 0; n < nclass; ++n) {
			last = nc[n] - 1;
			for (j = 0; j < sampsize[n]; ++j) {
			    ktmp = (int) (unif_rand() * (last+1));
			    k = class_idx[n][ktmp];
			    tmp = class_idx[n][last];
			    class_idx[n][last] = class_idx[n][ktmp];
			    class_idx[n][ktmp] = tmp;
			    last--;
			    tclasspop[cl[k] - 1] += wtt[k];
			    win[k] += wtt[k];
			    jin[k] = 1;
			}
		    }
		}
	    } else {  /* unstratified sampling */
		anyEmpty = 0;
		ntry = 0;
		do {
		    if (replace) {
			for (n = 0; n < *sampsize; ++n) {
			    k = unif_rand() * nsample;      
			    tclasspop[cl[k] - 1] += wtt[k];
			    win[k] += wtt[k];
			    jin[k] = 1;
			}
		    } else {
			for (n = 0; n < nsample; ++n) {
			    nind[n] = n;
			}
			last = nsample - 1;
			for (n = 0; n < *sampsize; ++n) {
			    ktmp = (int) (unif_rand() * (last+1));
			    k = nind[ktmp];
			    nind[ktmp] = nind[last];
			    nind[last] = k;
			    last--;
			    tclasspop[cl[k] - 1] += wtt[k];
			    win[k] += wtt[k];
			    jin[k] = 1;
			}
		    }
		    /* check if any class is missing in the sample */
		    for (n = 0; n < nclass; ++n) {
			if (tclasspop[n] == 0) anyEmpty = 1;
		    }
		    ntry++;
		} while (anyEmpty && ntry <= 5); 
	    }
      
	    /* Copy the original a matrix back. */ 
	    memcpy(a, at, sizeof(int) * mdim * nsample);
      	    modA(a, &nuse, nsample, mdim, cat, *maxcat, ncase, jin);
      
	    F77_CALL(buildtree)(a, b, cl, cat, &mdim, &nsample, &nclass, 
				treemap + 2*arrayindex, bestvar + arrayindex,
				bestsplit, bestsplitnext, tgini, 
				nodestatus + arrayindex, nodepop, 
				nodestart, classpop, tclasspop, tclasscat, 
				ta, nrnodes, idmove, &ndsize, ncase, parent, 
				jin, &mtry, varUsed, nodeclass + arrayindex, 
				ndbigtree + jb, win, wr, wc, wl, &mdim, 
				&nuse, mind);
	    /* if the "tree" has only the root node, start over */
	} while (ndbigtree[jb] == 1);
    
	Xtranslate(x, mdim, *nrnodes, nsample, bestvar + arrayindex, 
		   bestsplit, bestsplitnext, xbestsplit + arrayindex,
		   nodestatus + arrayindex, cat, ndbigtree[jb]);
    
	/*  Get test set error */
	if (*testdat) {
            predictClassTree(xts, ntest, mdim, zeroes, treemap + 2*arrayindex,
                             nodestatus + arrayindex, xbestsplit + arrayindex,
                             cbestsplit, bestvar + arrayindex, 
                             nodeclass + arrayindex, ndbigtree[jb], 
                             cat, nclass, jts, nodexts, *maxcat);
	    TestSetError(countts, jts, clts, outclts, ntest, nclass, jb+1,
			 errts + jb*(nclass+1), pid, *labelts, nclts, cut);
	}
    
	/*  Get out-of-bag predictions and errors. */
        predictClassTree(x, nsample, mdim, jin, treemap + 2*arrayindex,
                         nodestatus + arrayindex, xbestsplit + arrayindex,
                         cbestsplit, bestvar + arrayindex, 
                         nodeclass + arrayindex, ndbigtree[jb], 
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
#ifdef win32
	    R_FlushConsole();
	    R_ProcessEvents();
#endif
	    R_CheckUserInterrupt();
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
                    predictClassTree(x, nsample, mdim, jin, 
                                     treemap + 2*arrayindex,
                                     nodestatus + arrayindex, 
                                     xbestsplit + arrayindex,
                                     cbestsplit, 
                                     bestvar + arrayindex, 
                                     nodeclass + arrayindex,
                                     ndbigtree[jb], cat, nclass,
                                     jvr, nodex, *maxcat);

		    /* Count how often correct predictions are made with 
		       the modified data. */
		    for (n = 0; n < nsample; n++) {
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
			/* Restore the original data for that variable. */
		        x[m + n*mdim] = tx[n];
		    }
		    /* Accumulate decrease in proportions of correct 
		       predictions. */
		    for (n = 0; n < nclass; ++n) {
			if (nout[n] > 0) {
			    imprt[m + n*mdim] += 
				((double) (nright[n] - nrightimp[n])) / 
			        nout[n];
			    impsd[m + n*mdim] += 
				((double) (nright[n] - nrightimp[n]) * 
				 (nright[n] - nrightimp[n])) / nout[n]; 
			}
		    }
		    if (noutall > 0) {
			imprt[m + nclass*mdim] += 
				((double)(nrightall - nrightimpall)) / noutall;
			impsd[m + nclass*mdim] += 
				((double) (nrightall - nrightimpall) * 
				 (nrightall - nrightimpall)) / noutall; 
		    }
		}
	    }
	}

	/*  DO PROXIMITIES */
	if (iprox) {
	    for (n = 0; n < near; ++n) {
		for (k = n+1; k < near; ++k) {
		    if (oobprox) {
			if (jin[k] == 0 && jin[n] == 0) {
			    oobpair[k*near + n] ++;
			    oobpair[n*near + k] ++;
			    if (nodex[k] == nodex[n]) { 
				prox[k*near + n] += 1.0;
				prox[n*near + k] += 1.0;
			    }
			}
		    } else {
			if (nodex[k] == nodex[n]) {
			    prox[k*near + n] += 1.0;
			    prox[n*near + k] += 1.0;
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
		    for (k = 0; k < near; ++k) {
			if (nodexts[n] == nodex[k]) 
			    proxts[n + ntest * (k+ntest)] += 1.0; 
		    }
		}
	    } 
	}
	R_CheckUserInterrupt();
#ifdef win32
	R_ProcessEvents();
#endif
    }
    PutRNGstate();
      
    /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/
    if (iprox) {
	for (n = 0; n < near; ++n) {
	    for (k = n + 1; k < near; ++k) {
		if (oobprox) {
		    /* Only do the division if there are more than one 
		       instance!! */
		    if (oobpair[near*k + n] > 0) {
			prox[near*k + n] /= oobpair[near*k + n];
		    }
		} else {
		    prox[near*k + n] /= jbt;
		}
		prox[near*n + k] = prox[near*k + n];
	    }
	    prox[near*n + n] = 1.0;
	}
	if (*testdat) {
	    for (n = 0; n < ntest; ++n) 
		for (k = 0; k < ntest + nsample; ++k) 
		    proxts[ntest*k + n] /= jbt;
	}
    }
  
    /*  IMP DATA ++++++++++++++++++++++++++++++++++*/
    for (m = 0; m < mdim; m++) {
	tgini[m] = tgini[m] / jbt;
    }

    if (imp) {
	for (m = 0; m < mdim; ++m) {
	    /* casewise measures */
	    if (localImp) {
		for (n = 0; n < nsample; ++n) {
		    impmat[m + n*mdim] /= out[n];
		}
	    }
	    /* class-specific measures */
	    for (k = 0; k < nclass; ++k) {
	        av = imprt[m + k*mdim] / *ntree;
		impsd[m + k*mdim] = 
			sqrt(((impsd[m + k*mdim] / *ntree) - av*av) / *ntree);
		imprt[m + k*mdim] = av;
		/* imprt[m + k*mdim] = (se <= 0.0) ? -1000.0 - av : av / se; */
	    }
	    /* overall measures */
	    av = imprt[m + nclass*mdim] / *ntree;
	    impsd[m + nclass*mdim] = 
		    sqrt(((impsd[m + nclass*mdim] / *ntree) - av*av) / *ntree);
	    imprt[m + nclass*mdim] = av;
	    /* imprt[m + nclass*mdim] = (se <= 0.0) ? -1000.0-av : av/se; */
	    imprt[m + (nclass+1)*mdim] = tgini[m];
	}
    } else {
	for (m = 0; m < mdim; ++m) {
	    imprt[m] = tgini[m];
	}
    }
}

void runforest(int *mdim, int *ntest, int *nclass, int *maxcat,
	       int *nrnodes, int *jbt, 
	       double *xts, double *xbestsplit, double *pid, 
	       double *cutoff, double *countts, int *treemap,
	       int *nodestatus, int *cat, int *cbestsplit,
	       int *nodeclass, int *jts, int *jet, int *bestvar,
	       int *nodexts, int *ndbigtree, int *keepPred,
	       int *prox, double *proxmatrix, int *nodes) {
    int j, jb, n, n1, n2, idxNodes, idxSample1, idxSample2, *zeroes;
    double crit, cmax;

    zeroes = (int *) S_alloc(*ntest, sizeof(int));
    zeroInt(zeroes, *ntest);
    zeroDouble(countts, *nclass * *ntest);
    idxNodes = 0;
    idxSample1 = 0;
    idxSample2 = 0;

    for (jb = 0; jb < *jbt; ++jb) {
	/* predict by the jb-th tree */
        predictClassTree(xts, *ntest, *mdim, zeroes, treemap + 2*idxNodes,
			 nodestatus + idxNodes, xbestsplit + idxNodes,
			 cbestsplit, bestvar + idxNodes, nodeclass + idxNodes,
			 ndbigtree[jb], cat, *nclass,
			 jts + idxSample1, nodexts + idxSample2, *maxcat);

	/* accumulate votes: */
	for (n = 0; n < *ntest; ++n) {
	    countts[jts[n + idxSample1] - 1 + n * *nclass] += 1.0;
	}

	/* if desired, do proximities for this round */
	if (*prox) {
	    for (n1 = 0; n1 < *ntest; ++n1) {
		for (n2 = n1 + 1; n2 < *ntest; ++n2) {
		    if (nodexts[n1 + idxSample2] == nodexts[n2 + idxSample2]) {
			proxmatrix[n1 + n2 * *ntest] += 1.0;
		    }
		}
	    }
	}
	idxNodes += *nrnodes;
	if (*keepPred) idxSample1 += *ntest;
	if (*nodes)    idxSample2 += *ntest;
    }

    /* Aggregated prediction is the class with the maximum votes/cutoff */
    for (n = 0; n < *ntest; ++n) {
	cmax = 0.0;
	for (j = 0; j < *nclass; ++j) {
	    crit = (countts[j + n * *nclass] / *jbt) / cutoff[j];
	    if (crit > cmax) {
		jet[n] = j + 1;
		cmax = crit;
	    }
	    /* Break ties at random: */
	    if (crit == cmax) {
		if (unif_rand() > 0.5) {
		    jet[n] = j + 1;
		}
	    }
	}
    }

    /* if proximities requested, do the final adjustment 
       (division by number of trees) */
    if (*prox) {
	for (n1 = 0; n1 < *ntest; ++n1) {
	    for (n2 = n1 + 1; n2 < *ntest; ++n2) {
		proxmatrix[n1 + n2 * *ntest] /= *jbt;
		proxmatrix[n2 + n1 * *ntest] = proxmatrix[n1 + n2 * *ntest];
	    }
	    proxmatrix[n1 + n1 * *ntest] = 1.0;
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
    int j, n, noob, *noobcl;
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
	    for (j = 0; j < nclass; ++j) {
		qq = (((double) counttr[j + n*nclass]) / out[n]) / cutoff[j];
		if (j+1 != cl[n]) {
		    smax = (qq > smax) ? qq : smax;
		}
		/* if vote / cutoff is larger than current max, re-set max and 
		   change predicted class to the current class */
		if (qq > smaxtr) {
		    smaxtr = qq;
		    jest[n] = j+1;
		}
		/* break tie at random */
		if (qq == smaxtr) {
		    if (unif_rand() > 0.5) {
			smaxtr = qq;
			jest[n] = j+1;
		    }
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
    for (n = 1; n <= nclass; ++n) {
	errtr[n] /= noobcl[n-1];
    }
}


void TestSetError(double *countts, int *jts, int *clts, int *jet, int ntest,
		  int nclass, int nvote, double *errts, double *pid, 
		  int labelts, int *nclts, double *cutoff) {
    int j, n;
    double cmax, crit;
  
    for (n = 0; n < ntest; ++n) {
	countts[jts[n]-1 + n*nclass] += 1.0;
    }

    /*  Prediction is the class with the maximum votes */
    for (n = 0; n < ntest; ++n) {
	cmax=0.0;
	for (j = 0; j < nclass; ++j) {
	    crit = (countts[j + n*nclass] / nvote) / cutoff[j];
	    if (crit > cmax) {
		jet[n] = j+1;
		cmax = crit;
	    }

	    /*  Break ties at random: */
	    if (crit == cmax && unif_rand() > 0.5) {
		jet[n] = j+1;
		cmax = crit;
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
