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
#include <R.h>
#include "rf.h"

/*  Define the R RNG for use from Fortran. */
void F77_SUB(rrand)(double *r) { *r = unif_rand(); }

void rf(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat, 
	int *sampsize, int *Options, int *ntree, int *nvar,
	int *ipi, double *pi, double *cut, int *nodesize, 
        double *outlier, int *outcl, int *counttr, double *prox, 
	double *imprt, int *nrnodes, int *ndbigtree, 
	int *nodestatus, int *bestvar, int *treemap, int *nodeclass,
	double *xbestsplit, double *pid, double *errtr, 
	int *testdat, double *xts, int *clts, int *nts, double *countts,
	int *outclts, int *labelts, double *proxts, double *errts)
{
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
                  add a second class (for unsupervised RF)?
                     1: sampling from product of marginals
                     2: sampling from product of uniforms
                  assess variable importance?
                  calculate proximity?
                  calculate proximity based on OOB predictions?
                  calculate outlying measure?
                  how often to print output?
                  keep the forest for future prediction?
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
   *  prox:     matrix of proximity (if iprox=1)
   *  outlier:  measure of outlyingness (if noutlier=1)
   ******************************************************************/

  int nsample0, mdim, nclass, iaddcl, jbt, mtry, ntest, nsample, ndsize,
    mimp, nimp, near, nuse, noutall, nrightall, nrightimpall;
  int jb, j, n, m, mr, k, kt, itwo, arrayindex, imp, iprox, 
    oobprox, noutlier, keepf, trace, *nright, *nrightimp, *nout;

  int *out, *bestsplitnext, *bestsplit,
    *nodepop, *parent, *jin, *ndble, *nodex,
    *nodexts, *nodestart, *ta, *ncase, *jerr, *iv, *isort, *ncp, *clp,
    *jtr, *nc, *msum, *idmove, *jvr, /* *countimp, */
    *at, *a, *b, *cbestsplit, *mind, *jts, *oobpair;
  int **class_idx, *class_size, nboot;
  
  double errc, av=0, se=0;
  /* double maxgini = 0.0; */

  double *tgini, *v, *tx, *wl, *classpop, /* *errimp, *rimpmarg, *rmargin, */
    *tclasscat, *tclasspop, *win, *tp,
    *wc, *wr, *wtt, /* *diffmarg, *cntmarg, *rmissimp, */
    *tout, *tdx, *sm, *p, *iw, *sqsd, *sqsdall;

  iaddcl   = Options[0];
  imp      = Options[1];
  iprox    = Options[2];
  oobprox  = Options[3];
  noutlier = Options[4];
  trace    = Options[5]; 
  keepf    = Options[6];
  mdim     = dimx[0];
  nsample0 = dimx[1];
  nclass   = (*ncl==1) ? 2 : *ncl;
  nboot = 0;
  for (n = 0; n < nclass; ++n) {
    nboot += sampsize[n];
  }
  ndsize   = *nodesize;
  jbt      = *ntree;
  mtry     = *nvar;
  ntest    = *nts; 
  nsample = (iaddcl > 0) ? (nsample0 + nsample0) : nsample0;
  mimp = (imp == 1) ? mdim : 1;
  nimp = (imp == 1) ? nsample : 1;
  near = (iprox == 1) ? nsample0 : 1;
  if (trace == 0) trace = jbt + 1;


  tgini =      (double *) S_alloc(mdim, sizeof(double));
  v =          (double *) S_alloc(nsample, sizeof(double));
  tx =         (double *) S_alloc(nsample, sizeof(double));
  wl =         (double *) S_alloc(nclass, sizeof(double));
  classpop =   (double *) S_alloc(nclass* *nrnodes, sizeof(double));
  /* errimp =     (double *) S_alloc(mimp, sizeof(double)); */
  /* rimpmarg =   (double *) S_alloc(mdim*nsample, sizeof(double)); */
  tclasscat =  (double *) S_alloc(nclass*32, sizeof(double));
  tclasspop =  (double *) S_alloc(nclass, sizeof(double));
  /* rmargin =    (double *) S_alloc(nsample, sizeof(double)); */
  win =        (double *) S_alloc(nsample, sizeof(double));
  tp =         (double *) S_alloc(nsample, sizeof(double));
  wc =         (double *) S_alloc(nclass, sizeof(double));
  wr =         (double *) S_alloc(nclass, sizeof(double));
  wtt =        (double *) S_alloc(nsample, sizeof(double));
  /* diffmarg =   (double *) S_alloc(mdim, sizeof(double)); */
  /* cntmarg =    (double *) S_alloc(mdim, sizeof(double)); */
  /* rmissimp =   (double *) S_alloc(mimp, sizeof(double)); */
  tout =       (double *) S_alloc(near, sizeof(double));
  tdx =        (double *) S_alloc(nsample0, sizeof(double));
  sm =         (double *) S_alloc(nsample0, sizeof(double));
  p =          (double *) S_alloc(nsample0, sizeof(double));
  /* q =          (double *) S_alloc(nclass*nsample, sizeof(double)); */
  iw =         (double *) S_alloc(nsample, sizeof(double));
  sqsd =       (double *) S_alloc(mdim * nclass, sizeof(double));
  sqsdall =    (double *) S_alloc(mdim, sizeof(double));
  
  out =           (int *) S_alloc(nsample, sizeof(int));
  /* countimp =      (int *) S_alloc(nimp*mimp, sizeof(int)); */
  bestsplitnext = (int *) S_alloc(*nrnodes, sizeof(int));
  bestsplit =     (int *) S_alloc(*nrnodes, sizeof(int));
  nodepop =       (int *) S_alloc(*nrnodes, sizeof(int));
  parent =        (int *) S_alloc(*nrnodes, sizeof(int));
  jin =           (int *) S_alloc(nsample, sizeof(int));
  ndble =         (int *) S_alloc(nsample0, sizeof(int));
  nodex =         (int *) S_alloc(nsample, sizeof(int));
  nodexts =       (int *) S_alloc(ntest, sizeof(int));
  nodestart =     (int *) S_alloc(*nrnodes, sizeof(int));
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
  cbestsplit =    (int *) S_alloc(*maxcat * *nrnodes, sizeof(int));
  mind =          (int *) S_alloc(mdim, sizeof(int));
  nright =        (int *) S_alloc(nclass, sizeof(int));
  nrightimp =     (int *) S_alloc(nclass, sizeof(int));
  nout =          (int *) S_alloc(nclass, sizeof(int));
  if (oobprox) oobpair = (int *) S_alloc(near*near, sizeof(int));
  class_idx =     (int **) S_alloc(nclass, sizeof(int *));
  class_size =    (int *) S_alloc(nclass, sizeof(int));

  /* SET UP DATA TO ADD A CLASS++++++++++++++++++++++++++++++ */

  if(iaddcl >= 1) 
    F77_CALL(createclass)(x, cl, &nsample0, &nsample, &mdim, tdx, p,
                          sm, ndble, &iaddcl);
    
  /*    INITIALIZE FOR RUN */
  if(*testdat) {
    for (n = 0; n < nclass; ++n) {
      for (k = 0; k < ntest; ++k) countts[n + k*nclass] = 0.0;
    }
  }

  F77_CALL(zerm)(counttr, &nclass, &nsample);
  F77_CALL(zerv)(out, &nsample);
  F77_CALL(zervr)(tgini, &mdim);
  F77_CALL(zerv)(msum, &mdim);

  for(n = 0; n < jbt; ++n) errtr[n] = 0.0;

  if (*labelts) {
    for (n = 0; n < jbt; ++n) errts[n] = 0.0;
  }

  if (imp) {
    for (k = 0; k < nclass; ++k) {
      for (m = 0; m < mdim; ++m) {
	imprt[m + k*mdim] = 0.0;
	sqsd[m + k*mdim] = 0.0;
      }
    }
    for (n = 0; n < mdim; ++n) {
      sqsdall[n] = 0.0;
    }
  }

  if (iprox) {
    for (n = 0; n < nsample0; ++n)
      for (k = 0; k < nsample0; ++k)
         prox[k + n * nsample0] = 0.0;
    if (*testdat) {
      for (n = 0; n < ntest; ++n)
        for (k = 0; k < ntest + nsample0; ++k)
          proxts[n + k * ntest] = 0.0;
    }
  }
  
  F77_CALL(prep)(cl, &nsample, &nclass, ipi, pi, pid, nc, wtt);
  F77_CALL(makea)(x, &mdim, &nsample, cat, isort, v, at, b, &mdim);

  /* Create the array of pointers, each pointing to a vector of indices of 
     where data of each class is. */
  if (nboot) {
    for (n = 0; n < nclass; ++n) {
      class_idx[n] = (int *) S_alloc(nc[n], sizeof(int));
      class_size[n] = 0;
    }
    for (n = 0; n < nsample; ++n) {
      class_size[cl[n] - 1] ++;
      class_idx[cl[n] - 1][class_size[cl[n] - 1] - 1] = n;
    }
  }

  /*   START RUN   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  itwo = 2;
  GetRNGstate(); 

  for(jb = 0; jb < jbt; jb++) {

    arrayindex = (keepf == 1) ? jb *  *nrnodes : 0;

    F77_CALL(zerv)(nodestatus + arrayindex, nrnodes);
    F77_CALL(zerm)(treemap + 2*arrayindex, &itwo, nrnodes);
    F77_CALL(zervr)(xbestsplit + arrayindex, nrnodes);
    F77_CALL(zerv)(nodeclass + arrayindex, nrnodes);
    F77_CALL(zerv)(jin, &nsample);
    F77_CALL(zervr)(tclasspop, &nclass);
    F77_CALL(zervr)(win, &nsample);

    if (nboot) {
      for (n = 0; n < nclass; ++n) {
	for (j = 0; j < sampsize[n]; ++j) {
	  k = unif_rand() * nc[n];
	  k = class_idx[n][k];
	  tclasspop[cl[k] - 1] = tclasspop[cl[k] - 1] + wtt[k];
	  win[k] = win[k] + wtt[k];
	  jin[k] = 1;
	}
      }
    } else {
      for (n = 0; n < nsample; n++) {
	k = unif_rand() * nsample;      
	tclasspop[cl[k] - 1] = tclasspop[cl[k] - 1] + wtt[k];
	win[k] = win[k] + wtt[k];
	jin[k] = 1;
      }
    }
    
    /* Copy the original a matrix back. */
    for (n = 0; n < mdim*nsample; ++n) {
      a[n] = at[n];
    }

    F77_CALL(moda)(a, &nuse, &nsample, &mdim, cat, maxcat, ncase, jin, ta);
    
        F77_CALL(buildtree)(a, b, cl, cat, &mdim, &nsample, &nclass, 
             treemap + 2*arrayindex, 
             bestvar + arrayindex, bestsplit, bestsplitnext, tgini,
             nodestatus + arrayindex, nodepop, nodestart, classpop,
             tclasspop, tclasscat, ta, nrnodes, idmove,
             &ndsize, ncase, parent, jin, &mtry, iv,
             nodeclass + arrayindex, ndbigtree + jb, win, wr, wc,
             wl, &mdim, &nuse, mind); 

    F77_CALL(xtranslate)(x, &mdim, nrnodes, &nsample, bestvar + arrayindex, 
             bestsplit, bestsplitnext, xbestsplit + arrayindex,
             nodestatus + arrayindex, cat, ndbigtree + jb); 

    /*  Get test set error */
    if (*testdat) {
      F77_CALL(testreebag)(xts, &ntest, &mdim, treemap + 2*arrayindex,
			   nodestatus + arrayindex, xbestsplit + arrayindex,
			   cbestsplit, bestvar + arrayindex, 
			   nodeclass + arrayindex, nrnodes, ndbigtree + jb, 
			   cat, &nclass, jts, nodexts, maxcat);
      F77_CALL(comptserr)(countts, jts, clts, outclts, &ntest, &nclass,
               errts + jb, pid, labelts, cut);
    }
    
    /*  GET OUT-OF-BAG ESTIMATES */
    F77_CALL(testreebag)(x, &nsample, &mdim, treemap + 2*arrayindex, 
             nodestatus + arrayindex, xbestsplit + arrayindex, 
             cbestsplit, bestvar + arrayindex, 
             nodeclass + arrayindex, nrnodes, ndbigtree + jb, 
             cat, &nclass, jtr, nodex, maxcat); 

    for (n = 0; n < nclass; ++n) {
      nout[n] = 0;
    }
    noutall = 0;
    for (n = 0; n < nsample; ++n) {
      if (jin[n] == 0) {
	/* increment the OOB votes */
        counttr[n*nclass + jtr[n] - 1] ++;
	/* count number of times a case is OOB */
        out[n]++;
	/* count number of OOB cases in the current iteration. */
	nout[cl[n] - 1]++;
	noutall++;
      }
    }

    F77_CALL(oob)(&nsample, &nclass, jin, cl, jtr, jerr, counttr, out,
             errtr + jb, &errc, outcl, wtt, cut);

    if ((jb+1) % trace == 0) {
      Rprintf("%4i: OOB error rate=%5.2f%%\t", jb+1, 100.0*errtr[jb]);
      if (*labelts) 
	Rprintf("Test set error rate=%5.2f%%", 100.0*errts[jb]);
      Rprintf("\n");
    }

    /*  DO VARIABLE IMPORTANCE  */
    if (imp) {
      nrightall = 0;
      F77_CALL(zerv)(iv, &mdim);
      for (kt = 0; kt < *(ndbigtree + jb); kt++) {
        if (nodestatus[kt + arrayindex] != -1) {
          iv[bestvar[kt + arrayindex] - 1] = 1;
          msum[bestvar[kt + arrayindex] - 1] ++;
	}
      }
      /* Count the number of correct prediction by the current tree
	 among the OOB samples, by class. */
      for (n = 0; n < nclass; ++n) {
	nright[n] = 0;
      }

      for (n = 0; n < nsample; ++n) {
	if (jin[n] == 0 && jtr[n] == cl[n]) {
	    nright[cl[n] - 1]++;
	    nrightall++;
	}
      }

      for (mr = 1; mr <= mdim; mr++) { 
        if (iv[mr-1] == 1) {
	  nrightimpall = 0;
	  /* Permute the mr-th variable. */
          F77_CALL(permobmr)(&mr, x, tp, tx, jin, &nsample, &mdim);
	  /* Predict the modified data using the current tree. */
          F77_CALL(testreebag)(x, &nsample, &mdim, treemap + 2*arrayindex, 
                               nodestatus + arrayindex,
                               xbestsplit + arrayindex, cbestsplit, 
                               bestvar + arrayindex, nodeclass + arrayindex, 
                               nrnodes, ndbigtree + jb, cat,
                               &nclass, jvr, nodex, maxcat); 
	  for (n = 0; n < nclass; ++n) {
	    nrightimp[n] = 0;
	  }
	  /* Count how often correct predictions are made with the modified
	     data. */
          for (n = 0; n < nsample; n++) {
            if (jin[n] == 0 && cl[n] == jvr[n]) {
	      nrightimp[cl[n] - 1]++;
	      nrightimpall++;
            }
          }
	  /* Accumulate decrease in proportions of correct predictions. */
	  for (n = 0; n < nclass; ++n) {
	    /*	    Rprintf("%d %d %d ", nright[n], nrightimp[n], nout[n]);*/
	    if (nout[n] > 0) {
	      imprt[(mr - 1) + n*mdim] += 
		((double) (nright[n] - nrightimp[n])) / nout[n];
	      sqsd[(mr - 1) + n*mdim]  += 
		((double) (nright[n] - nrightimp[n]) * 
		 (nright[n] - nrightimp[n])) / nout[n]; 
	    }
	  }
	  imprt[(mr-1) + nclass*mdim] += ((double)(nrightall - nrightimpall)) /
	    noutall;
	  sqsdall[mr-1] += ((double) (nrightall - nrightimpall) * 
		 (nrightall - nrightimpall)) / noutall; 

	  /* Put back the original data for that variable. */
          for (n = 0; n < nsample; n++) {
            x[mr-1 + n*mdim] = tx[n];
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
            if (nodexts[n] == nodex[k]) proxts[n + ntest * (k+ntest)] += 1.0; 
          } 
        }
      } 
    }
  }
  PutRNGstate();
      
  /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/  
  if (iprox) {
    for (n = 0; n < near; ++n) {
      for (k = n + 1; k < near; ++k) {
	if (oobprox) {
	  /* Only do the division if there are more than one instance!! */
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

  /*  OUTLIER DATA +++++++++++++++++++++++++++++*/
  if (noutlier) {
    F77_CALL(locateout)(prox, cl, &near, &nsample, &nclass, ncp, &iaddcl,
             outlier, tout, isort, clp);
  }
  
  /*  IMP DATA ++++++++++++++++++++++++++++++++++*/
  for (m = 0; m < mdim; m++) {
    tgini[m] = tgini[m] / jbt;
  }
  if (imp) {
    /*
    F77_CALL(finishimp)(rmissimp, countimp, out, cl, &nclass, &mdim,
                        &nsample, errimp, rimpmarg, diffmarg, cntmarg, 
                        rmargin, counttr, outcl, errtr + jbt - 1, cut);
    */
    for (m = 0; m < mdim; ++m) {
      /* class-specific measures */
      for (k = 0; k < nclass; ++k) {
	av = imprt[m + k*mdim] / *ntree;
	/* se = sqrt( ((sqsd[m + k*mdim] / *ntree) - av * av) / *ntree ); */
	imprt[m + k*mdim] = av; /* / se; */
      }
      /* overall measures */
      av = imprt[m + nclass*mdim] / *ntree;
      se = sqrt( ((sqsdall[m] / *ntree) - av * av) / *ntree );
      imprt[m + nclass*mdim] = av; /* / se; */
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
	       int *prox, double *proxmatrix)
{      
  int i, j, jb, n, n1, n2, idxNodes, idxSample;
  double crit, cmax;

  for (i = 0; i < *nclass * *ntest; ++i) {
    countts[i] = 0.0;
  }

  for (jb = 0; jb < *jbt; ++jb) {

    idxNodes = jb * *nrnodes;
    /* Offset for prediction by individual tree: if prediction is to be 
       kept, advance the index by the number of cases. */
    idxSample = *keepPred ? jb * *ntest : 0;

    /* predict by the jb-th tree */
    F77_CALL(testreebag)(xts, ntest, mdim, treemap + 2*idxNodes,
			 nodestatus + idxNodes, xbestsplit + idxNodes,
			 cbestsplit + idxNodes, bestvar + idxNodes, 
			 nodeclass + idxNodes, nrnodes, ndbigtree + jb, 
			 cat, nclass, jts + idxSample, nodexts, maxcat);

    /* accumulate votes: */
    for (n = 0; n < *ntest; ++n) {
      countts[jts[n + idxSample] - 1 + n * *nclass] += 1.0;
    }

    /* if desired, do proximities for this round */
    if (*prox) {
      for (n1 = 0; n1 < *ntest; ++n1) {
	for (n2 = n1 + 1; n2 < *ntest; ++n2) {
	  if (nodexts[n1] == nodexts[n2]) {
	    proxmatrix[n1 + n2 * *ntest] += 1.0;
	  }
	}
      }
    }
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
