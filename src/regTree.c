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

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's 
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose. 
 *      It comes with no guarantee.  
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"

void regTree(double *x, double *y, int mdim, int nsample, 
	     int *treemap, double *upper, double *avnode, int *nodestatus, 
	     int nrnodes, int nthsize, int mtry, int *mbest, int *cat, 
	     double *tgini) {
  /* int mdim, nsample, nrnodes, nthsize, mtry; */
    int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
    int ndstart, ndend, ndendl, nodecnt, jstat, msplit;
    double d, ss, av, decsplit, ubest, sumnode;

    /*    mtry = *mtr;
    mdim = *mdi;
    nsample = *nsampl;
    nrnodes = *nrnode;
    nthsize = *nthsiz;*/

    nodestart = (int *) R_alloc(nrnodes, sizeof(int));
    nodepop   = (int *) R_alloc(nrnodes, sizeof(int));
    
    /* initialize some arrays for the tree */
    for (i = 0; i < nrnodes; ++i) {
	nodestatus[i] = 0;
	nodestart[i] = 0;
	nodepop[i] = 0;
	avnode[i] = 0;
    }
    
    jdex = (int *) R_alloc(nsample, sizeof(int));
    for (i = 1; i <= nsample; ++i) {
	jdex[i-1] = i;
    }
    for (i=0; i < mdim; ++i) tgini[i] = 0.0;

    ncur = 0;
    nodestart[0] = 0;
    nodepop[0] = nsample;
    nodestatus[0] = 2;
    
    /* compute mean and sum of squares for Y */
    av=0.0;
    ss=0.0;
    for (i=0; i < nsample; ++i) {
	d = y[jdex[i] - 1];
	ss += i * (av - d) * (av - d) / (i + 1);
	av = (i * av + d) / (i+1);
    }
    avnode[0] = av;
    
    /* start main loop */
    for (k = 0; k < nrnodes; ++k) {
	if (k > ncur) break;
	/* Rprintf("k=%i, status=%i, ncur=%i\n", k, nodestatus[k], ncur); */
	if (nodestatus[k] != 2) continue;
	
	/* initialize for next call to findbestsplit */         
	ndstart = nodestart[k];
	ndend = ndstart + nodepop[k] - 1;
	nodecnt = nodepop[k]; 
	sumnode = nodecnt * avnode[k];
	jstat = 0;
	decsplit = 0.0;
	
	findBestSplit(x, jdex, y, mdim, nsample,
		      ndstart, ndend, &msplit, &decsplit, &ubest, 
		      &ndendl, &jstat, mtry, sumnode, nodecnt, cat);
	if (jstat == 1) {
	    /* Node is terminal: Mark it as such and move on to the next. */
	    nodestatus[k] = -1;
	    continue;
	} else {
	    /* Found the best split. */
	    mbest[k] = msplit;
	    upper[k] = ubest;
/*	    bestcrit[k] = decsplit; */
	}
	tgini[msplit-1] += decsplit;
	
	/* leftnode no.= ncur+1, rightnode no. = ncur+2. */
	nodepop[ncur + 1] = ndendl - ndstart + 1;
	nodepop[ncur + 2] = ndend - ndendl;
	nodestart[ncur + 1] = ndstart;
	nodestart[ncur + 2] = ndendl + 1;
    
	/* compute mean and sum of squares for the left daughter node */
	av = 0.0;
	ss = 0.0;
	for (j = ndstart; j <= ndendl; ++j) {
	    d = y[jdex[j]-1];
	    m = j - ndstart;
	    ss += m * (av - d) * (av - d) / (m + 1);
	    av = (m * av + d) / (m+1);
	}
 	/* Rprintf("node %i, av=%f, count=%i, start=%i, end=%i\n", ncur+1, av, 
	   nodepop[ncur+1], ndstart, ndendl); 
	   for (j = ndstart; j <= ndendl; ++j) {
	   Rprintf("j=%i, y[%i]=%f\n", j, jdex[j] - 1, y[jdex[j]-1]);
	  } */
	avnode[ncur+1] = av;
	nodestatus[ncur+1] = 2;
	if (nodepop[ncur + 1] <= nthsize) {
	    nodestatus[ncur + 1] = -1;
	}
	
	/* compute mean and sum of squares for the right daughter node */
	av = 0.0;
	ss = 0.0;
	for (j = ndendl + 1; j <= ndend; ++j) {
	    d = y[jdex[j]-1];
	    m = j - (ndendl + 1);
	    ss += m * (av - d) * (av - d) / (m + 1);
	    av = (m * av + d) / (m + 1);
	}
 	/* Rprintf("node %i, av=%f, count=%i, start=%i, end=%i\n", ncur+2, av, 
	   nodepop[ncur+2], ndendl+1, ndend); 
	   for (j = ndendl + 1; j <= ndend; ++j) {
	   Rprintf("y[%i]=%f\n", jdex[j] - 1, y[jdex[j]-1]);
	  } */
	avnode[ncur + 2] = av;
	nodestatus[ncur+2] = 2;
	if (nodepop[ncur + 2] <= nthsize) {
	    nodestatus[ncur + 2] = -1;
	}
	
	/* map the daughter nodes */
	treemap[k * 2] = ncur + 1 + 1;
	treemap[1 + k * 2] = ncur + 2 + 1;
	nodestatus[k] = 1;
	/* Augment the tree by two nodes. */
	ncur += 2;
	/* if tree size limit has been reach, stop */
	if (ncur >= nrnodes) break;
    }
/*
    Free(jdex);
    Free(nodepop);
    Free(nodestart);
    Free(bestcrit);*/
}

/*--------------------------------------------------------------*/
void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample, 
		   int ndstart, int ndend, int *msplit, double *decsplit, 
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat) {
    int last, ncat[32], icat[32], non, lc, nl, nr, npopl, npopr, ic;
    int i, j, kv, l, *mind, *ncase;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], ubestt;
    double crit, critmax, critvar, suml, sumr, d;

    /*
    ut = (double *) Calloc(nsample, double);
    xt = (double *) Calloc(nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    */
    ut = (double *) S_alloc(nsample, sizeof(double));
    xt = (double *) S_alloc(nsample, sizeof(double));
    v  = (double *) S_alloc(nsample, sizeof(double));
    yl = (double *) S_alloc(nsample, sizeof(double));
    mind  = (int *) S_alloc(mdim, sizeof(int));
    ncase = (int *) S_alloc(nsample, sizeof(int));
    
    /* START BIG LOOP */
    *msplit = -1;
    *decsplit = 0.0;
    critmax = 0.0;
    ubestt = 0.0;
    non=0;
    for (i=0; i < mdim; ++i) {
	mind[i] = i;
    }
    
    last = mdim - 1;
    for (i = 0; i < mtry; ++i) {
	critvar = 0.0;
	j = (int) (unif_rand() * (last+1));
	kv = mind[j];
	mind[j] = mind[last];
	mind[last] = kv;
	last--;

	lc = cat[kv];
	if (lc == 1) {
	    /* numeric variable */
	    for (j = ndstart; j <= ndend; ++j) {
		 xt[j] = x[kv + (jdex[j]-1) * mdim];
		 yl[j] = y[jdex[j]-1];
		 /* Rprintf("j=%i, xt=%f, yl=%f\n", j, xt[j], yl[j]); */
	    } 
	} else {
	    /* categorical variable */
	    for (j=0; j < 32; ++j) {
		sumcat[j] = 0.0;
		ncat[j] = 0;
	    }
	    for (j = ndstart; j <= ndend; ++j) {
		l = (int) x[kv + (jdex[j]-1) * mdim];
		d = y[jdex[j]-1];
		sumcat[l-1] += d;
		ncat[l-1] ++;
	    }
	    for (j = 0; j < lc; ++j) {
		if (ncat[j] > 0) {
		    avcat[j] = sumcat[j] / ncat[j];
		} else {
		    avcat[j] = 0.0;
		}
	    }
	    for (j = 0; j < nsample; ++j) {
		xt[j] = avcat[(int) x[kv + (jdex[j]-1) * mdim]];
		yl[j] = y[jdex[j]-1];
	    }
	}
	for (j = ndstart; j <= ndend; ++j) {
	    v[j] = xt[j];
	}
	
	for (j = 1; j <= nsample; ++j) {
	    ncase[j-1] = j;
	}
	
	R_qsort_I(v, ncase, ndstart+1, ndend+1);
	for (j = ndstart; j <= ndend; ++j) {
	}
	if (v[ndstart] >= v[ndend]) {
	    non++;
	    if (non >= 3*mtry) {
		*jstat = 1;
		return;
	    }
	}
	    
	/* ncase(n)=case number of v nth from bottom */
	/* Start from the right and search to the left. */
	suml = 0.0;
	sumr = sumnode;
	npopl = 0;
	npopr = nodecnt;
	crit = 0.0;
	/* Rprintf("start=%i   end=%i\n", ndstart, ndend); */
	/* Search through the "gaps" in the x-variable. */
	for (j = ndstart; j <= ndend-1; ++j) {
	  /* Rprintf("yl[ncase[%i]]=%f\n", j, yl[ncase[j]-1]); */
	    d = yl[ncase[j]-1];
	    suml += d;
	    sumr -= d;
	    npopl++;
	    npopr--;
	    if (v[j] < v[j+1]) {
		crit = (suml * suml / npopl) + (sumr * sumr / npopr);
		if (crit > critvar) {
		    ubestt = (v[j] + v[j+1]) / 2.0;
		    critvar = crit;
		    /* nbestt = j; */
		}
	    }
	}
	if (critvar > critmax) {
	    *ubest = ubestt;
	    /* nbest = nbestt; */
	    *msplit = kv + 1;
	    critmax = critvar;
	    for (j = ndstart; j <= ndend; ++j) {
		ut[j] = xt[j];
	    }
	    if (cat[kv] > 1) {
		ic = cat[kv];
		for (j = 0; j < ic; ++j) {
		    tavcat[j] = avcat[j];
		}
	    }
	}
	/* Rprintf("var %i: crit=%f, critmax=%f\n", kv+1, critvar, critmax);*/
    }
    *decsplit = critmax - (sumnode * sumnode / nodecnt);
    
    /* If best split can not be found, set to terminal node and return. */
    if (*msplit == -1) {
	*jstat = 1;
	return;
    }

    nl = ndstart;
    for (j = ndstart; j <= ndend; ++j) {
	if (ut[j] <= *ubest) {
	    nl++;
	    ncase[nl-1] = jdex[j];
	}
    }
    *ndendl = imax2(nl - 1, ndstart);
    nr = *ndendl + 1;
    for (j = ndstart; j <= ndend; ++j) {
	if (ut[j] > *ubest) {
	    if (nr >= nsample) break;
	    nr++;
	    ncase[nr-1] = jdex[j];
	}
    } 
    if (*ndendl >= ndend) *ndendl = ndend - 1; 

    for (j = ndstart; j <= ndend; ++j) {
      jdex[j] = ncase[j];
    }

    lc = cat[*msplit-1];
    if (lc > 1) {
	for (j = 0; j < lc; ++j) {
	    icat[j] =  (tavcat[j] < *ubest) ? 1 : 0;
	}
	*ubest = pack(lc, icat);
    }
    /*
      Free(ncase);
      Free(mind);
      Free(v);
      Free(yl);
      Free(xt);
      Free(ut); 
    */
}

/*====================================================================*/
void predictRegTree(double *x, int nsample, int mdim, int *doPred,
		    int *treemap, int *nodestatus, int nrnodes, 
		    int ndbigtree, double *ypred, double *split, 
		    double *nodepred, int *bestvar, int *cat, int *nodex) {
    int icat[32], i, j, k, kt, m, mm, lc;

    zeroInt(nodex, nsample);

    for (i = 0; i < nsample; ++i) {
	if (doPred[i] > 0) continue; /* skip to the next case */
	kt = 0;
	for (j = 0; j < ndbigtree; ++j) {
	    if (nodestatus[kt] == -1) {
		/* terminal node: assign prediction and move on to next */
		ypred[i] = nodepred[kt];
		nodex[i] = kt + 1;
		break;
	    }
	    m = bestvar[kt] - 1;
	    lc = cat[m];
	    if (lc == 1) {
		kt = (x[m + i*mdim] <= split[kt]) ? 
		    treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    } else {
		/* Should change this part to using bit manipulations... */
		mm = (int) split[kt];
		unpack(lc, mm, icat);
		k = (int) x[m + i * mdim];
		kt = icat[k-1] ? treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    }
	}
    }
}


double pack(int l, int *icat) {
    /* icat is a binary integer with ones for categories going left 
     * and zeroes for those going right.  The sub returns npack- the integer */
    int k;
    double pack = 0.0;

    for (k = 0; k < l; ++k) {
	if (icat[k]) pack += R_pow_di(2.0, k);
    }
    return(pack);
}


void unpack(int l, int npack, int *icat) {
/*      
 * npack is a long integer.  The sub. returns icat, an integer of zeroes and
 * ones corresponding to the coefficients in the binary expansion of npack.
 */   
    int i;
    for (i = 0; i < 32; ++i) {
	icat[i] = 0;
    }
    icat[0] = npack % 2;
    for (i = 1; i < l; ++i) {
	npack = (npack - icat[i-1]) / 2;
	icat[i] = npack % 2;
    }
}
