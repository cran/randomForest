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
#include <Rmath.h>
#include "rf.h"

#ifdef C_CLASSTREE
void classTree(int *a, int *b, int *class, int *cat, int mdim, int nsample, 
               int nclass, int *treemap, int *bestvar, double *bestsplit,
               double *bestsplitnext, double *tgini, int *nodeStatus,
               int *nodePop, int *nodeStart, double *tclassPop, int maxNodes,
               int nodeSize, int *ncase, int *inBag, int mTry, int *varUsed,
               int *nodeClass, int *treeSize, double *win) {
/* Buildtree consists of repeated calls to two subroutines, Findbestsplit
   and Movedata.  Findbestsplit does just that--it finds the best split of
   the current node.  Movedata moves the data in the split node right and
   left so that the data corresponding to each child node is contiguous.
   The buildtree bookkeeping is different from that in Friedman's original
   CART program.  ncur is the total number of nodes to date.
   nodeStatus(k)=1 if the kth node has been split.  nodeStatus(k)=2 if the
   node exists but has not yet been split, and =-1 of the node is terminal.
   A node is terminal if its size is below a threshold value, or if it is
   all one class, or if all the x-values are equal.  If the current node k
   is split, then its children are numbered ncur+1 (left), and
   ncur+2(right), ncur increases to ncur+2 and the next node to be split is
   numbered k+1.  When no more nodes can be split, buildtree returns to the
   main program.
*/
/*
  integer a(mdim,nsample),cl(nsample),cat(mdim),
  treemap(2,numNodes),bestvar(numNodes),
          bestsplit(numNodes), nodeStatus(numNodes),ta(nsample),
          nodePop(numNodes),nodeStart(numNodes),
          bestsplitnext(numNodes),idmove(nsample),
          ncase(nsample),parent(numNodes),b(mdim,nsample),
          jin(nsample),iv(mred),nodeclass(numNodes),mind(mred)
      
      
      double precision tclasspop(nclass),classpop(nclass,numNodes),
     1     tclasscat(nclass,32),win(nsample),wr(nclass),wc(nclass),
     1     wl(nclass),tgini(mdim), xrand
 */          
    int msplit = 0, i, j;      
    zeroInt(nodeStatus, maxNodes);
    zeroInt(nodeStart, maxNodes);
    zeroInt(nodePop, maxNodes);
    zeroDouble(classPop, nclass * maxNodes);
      
    for (i = 0; i < nclass; ++i) classPop[i] = tclassPop[i];
    ncur = 1;
    nodeStart[0] = 1;
    nodePop[0] = *nuse;
    nodeStatus[0] = NODE_TOSPLIT;
    /* 2: not split yet, 1: split, -1: terminal */      
    /* start main loop */
    for (i = 0; i < numNodes; ++i) {
        if (i > ncur - 1) break;
        if (nodeStatus[i] != NODE_TOSPLIT) continue; 
        /* initialize for next call to findbestsplit */
        ndstart = nodeStart[i];
        ndend = ndstart + nodePop[i] - 1;
        for (j = 0; j < nclass; ++j) {
            tclassPop[j] = classPop[j + i * nclass];
        }
        jstat = 0;
        F77_CALL(findbestsplit)(a, b, cl, mdim, nsample, nclass, cat,
                                ndstart, ndend, tclassPop, tclasscat,
                                &msplit, &decsplit, &nbest, ncase, &jstat,
                                inBag, mTry, win, wr, wc, wl, mred, i, mind);
        if (jstat == 1) { 
            nodeStatus[i] = NODE_TERMINAL;
            continue;
        } else { 
            bestvar[i] = msplit;
            varUsed[msplit - 1] = 1;
            tgini[msplit - 1] += decsplit;
            if (cat[msplit-1] == 1) {
                bestsplit[i] = a[msplit - 1  + nbest * mdim];
                bestsplitnext[i] = a[msplit - 1 + (nbest + 1) * mdim];
            } else {
                  bestsplit[i] = nbest;
                  bestsplitnext[i] = 0;
            }
        }          
        F77_CALL(movedata)(a, ta, mdim, nsample, ndstart, ndend, idmove,
                           ncase, msplit, cat, nbest, ndendl);
        /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
        nodePop[ncur+1] = ndendl - ndstart + 1;
        nodePop[ncur+2] = ndend - ndendl;
        nodeStart[ncur+1] = ndstart;
        nodeStart[ncur+2] = ndendl + 1;
        /* find class populations in both nodes */
        for (n = ndstart; n <= ndendl; ++n) {
            nc = ncase[n];
            j = class[nc-1];
            classPop[j - 1 + (ncur+1)*mdim] += win[nc - 1];
        }
        for (n = ndendl + 1; n <= ndend; ++n) {
            nc = ncase[n];
            j = cl[nc - 1];
            classPop[j - 1 + (ncur+2) * mdim] += win[nc - 1];
        }
        /* check on nodeStatus */
        nodeStatus[ncur + 1] = NODE_TOSPLIT;
        nodeStatus[ncur + 2] = NODE_TOSPLIT;
        if (nodePop[ncur + 1] <= ndsize) nodeStatus[ncur+1] = NODE_TERMINAL;
        if (nodePop[ncur + 2] <= ndsize) nodeStatus[ncur+2] = NODE_TERMINAL;
        popt1 = 0;
        popt2 = 0;
        for (j = 0; j < nclass; ++j) {
            popt1 += classPop[j + (ncur+1) * mdim];
            popt2 += classPop[j + (ncur+2) * mdim];
        }
        for (j = 0; j < nclass; ++j) {
            if (classPop[j + (ncur+1) * mdim] == popt1) 
                nodeStatus[ncur+1] = NODE_TERMINAL;
            if (classPop[j + (ncur+2) * mdim] == popt2) 
                nodeStatus[ncur+2] = NODE_TERMINAL;
        }
        
        treemap[i * 2] = ncur + 1;
        treemap[1 + i * 2] = ncur + 2;
        nodeStatus[i] = NODE_INTERIOR;
        ncur += 2;
        if (ncur >= numNodes) break;
    }
    ndbigtree = numNodes;
    for (k = numNodes-1; k >= 0; --k) {
        if (nodeStatus[k] == 0) ndbigtree--;
        if (nodeStatus[k] == NODE_TOSPLIT) nodeStatus[k] = NODE_TERMINAL;
    }
    for (k = 0; k < ndbigtree; ++k) {
        if (nodeStatus[k] == NODE_TERMINAL) {
            pp = 0;
            for (j = 0; j < nclass; ++j) {
                if (classPop[j + k * nclass] > pp) {
                    nodeClass[k] = j;
                    pp = classPop[j + k * nclass];
                }
                /* Break ties at random: */
                if (classPop[j + k * nclass] == pp && unif_rand() > 0.5) {
                    nodeClass[k] = j;
                    pp = classPop[j + k * nclass];
                }
            }
        }
    }
}
#endif /* C_CLASSTREE */

void predictClassTree(double *x, int n, int mdim, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *bestvar, int *nodeclass,
		      int treeSize, int *cat, int nclass,
		      int *jts, int *nodex, int maxcat) {
    int m, npack, i, j, k, *cbestsplit;

    /* decode the categorical splits */
    if (maxcat > 1) {
        cbestsplit = (int *) Calloc(maxcat * treeSize, int);
        zeroInt(cbestsplit, maxcat * treeSize);
        for (i = 0; i < treeSize; ++i) {
            if (nodestatus[i] != NODE_TERMINAL) {
                if (cat[bestvar[i] - 1] > 1) {
                    npack = (int) xbestsplit[i];
                    /* unpack `npack' into bits */
                    for (j = 0; npack; npack >>= 1, ++j) {
                        cbestsplit[j + i*maxcat] = npack & 01;
                    }
                }
            }
        }
    }
    for (i = 0; i < n; ++i) {
	k = 0;
	while (nodestatus[k] != NODE_TERMINAL) {
            m = bestvar[k] - 1;
            if (cat[m] == 1) {
  	        /* Split by a numerical predictor */
	        k = (x[m + i * mdim] <= xbestsplit[k]) ?
		    treemap[k * 2] - 1 : treemap[1 + k * 2] - 1;
	    } else {
	        /* Split by a categorical predictor */
	        k = cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] ?
		    treemap[k * 2] - 1 : treemap[1 + k * 2] - 1;
	    }
	}
	/* Terminal node: assign class label */
	jts[i] = nodeclass[k];
	nodex[i] = k + 1;
    }
    if (maxcat > 1) Free(cbestsplit);
}

void F77_NAME(catmax)(double *pdo, double *tclasscat, double *tclasspop,
                       int *nclass, int *lcat, int *ncatsp, 
                       double *critmax, int *nhit, int *maxcat, int *ncmax,
                       int *ncsplit) {
/* This finds the best split of a categorical variable with lcat 
   categories and nclass classes, where tclasscat(j, k) is the number 
   of cases in class j with category value k. The method uses an 
   exhaustive search over all partitions of the category values if the 
   number of categories is 10 or fewer.  Otherwise ncsplit randomly 
   selected splits are tested and best used. */
    int j, k, n, icat[32], nsplit;
    double pln, pld, prn, tdec, *tmpclass;
    
    tmpclass = (double *) Calloc(*nclass, double);
    *nhit = 0;
    nsplit = *lcat > *ncmax ? 
        *ncsplit : (int) pow(2.0, (double) *lcat - 1) - 1;

    for (n = 0; n < nsplit; ++n) {
        zeroInt(icat, 32);
        if (*lcat > *ncmax) {
            /* Generate random split.
               TODO: consider changing to generating random bits with more
               efficient algorithm */
            for (j = 0; j < *lcat; ++j) icat[j] = unif_rand() > 0.5 ? 1 : 0;
        } else {
            unpack(n, icat);
        }
        for (j = 0; j < *nclass; ++j) {
            tmpclass[j] = 0;
            for (k = 1; k < *lcat; ++k) {
                if (icat[k]) tmpclass[j] += tclasscat[j + k * *nclass];
            }
        }
        pln = 0.0;
        pld = 0.0;
        for (j = 0; j < *nclass; ++j) {
            pln += tmpclass[j] * tmpclass[j];
            pld += tmpclass[j];
        }
        prn = 0.0;
        for (j = 0; j < *nclass; ++j) {
            tmpclass[j] = tclasspop[j] - tmpclass[j];
            prn += tmpclass[j] * tmpclass[j];
        }
        tdec = (pln / pld) + (prn / (*pdo - pld));
        if (tdec > *critmax) {
            *critmax = tdec;
            *nhit = 1;
            *ncatsp = *lcat > *ncmax ? pack(*lcat, icat) : n;
        }
    }
    Free(tmpclass);
}

/* Find best split of with categorical variable when there are two classes */
void F77_NAME(catmaxb)(double *pdo, double *tclasscat, double *tclasspop,
                       int *nclass, int *lcat, int *nbest, double *critmax, 
                       int *nhit, int *maxcat, double *dn) {

    double xc[32], cp[32], cm[32];
    int kcat[32];
    int i, j;
    double bestsplit=0.0, rrd, rld, rln, rrn, crit;

    *nhit = 0;
    for (i = 0; i < *lcat; ++i) {
        xc[i] = dn[i] ? tclasscat[i * *nclass] / dn[i] : 0;
        kcat[i] = i + 1;
    }
    R_qsort_I(xc, kcat, 1, *lcat);
    for (i = 0; i < *nclass; ++i) {
        cp[i] = 0;
        cm[i] = tclasspop[i];
    }
    rrd = *pdo;
    rld = 0.0;
    for (i = 0; i < *lcat - 1; ++i) {
        rld += dn[kcat[i]];
        rrd -= dn[kcat[i]];
        rln = 0.0;
        rrn = 0.0;
        for (j = 0; j < *nclass; ++j) {
            cp[j] += tclasscat[j + kcat[i] * *nclass];
            cm[j] -= tclasscat[j + kcat[i] * *nclass];
            rln += cp[j] * cp[j];
            rrn -= cm[j] * cm[j];
        }
        if (xc[i] < xc[i + 1]) {
            if (fmin2(rrd, rld) > 1.0) {
                crit = (rln / rld) + (rrn / rrd);
                if (crit > *critmax) {
                    *critmax = crit;
                    bestsplit = .5 * (xc[i] + xc[i + 1]);
                    *nhit = 1;
                }
            }
        }
    }
    if (*nhit == 1) {
        zeroInt(kcat, *maxcat);
        for (i = 0; i < *lcat; ++i) {
            xc[i] = dn[i] ? tclasscat[i * *nclass] / dn[i] : 0.0;
            kcat[i] = xc[i] < bestsplit ? 1 : 0;
        }
        *nbest = pack(*lcat, kcat);
    }
}
