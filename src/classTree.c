#include <R.h>
#include "rf.h"
/*
void classTree(int *a, int *b, int *class, int *cat, int mdim, int nsample, 
               int nclass, int *treemap, int *bestvar, double *bestsplit,
               double *bestsplitnext, double *tgini, int *nodeStatus,
               int *nodePop, int *nodeStart, double *tclassPop, int numNodes,
               int nodeSize, int *ncase, int *inBag, int mTry, int *varUsed,
               int *nodeClass, int *totNodes, double *win) {

  c Buildtree consists of repeated calls to two subroutines, Findbestsplit
  c and Movedata.  Findbestsplit does just that--it finds the best split of
  c the current node.  Movedata moves the data in the split node right and
  c left so that the data corresponding to each child node is contiguous.
  c The buildtree bookkeeping is different from that in Friedman's original
  c CART program.  ncur is the total number of nodes to date.
  c nodestatus(k)=1 if the kth node has been split.  nodestatus(k)=2 if the
  c node exists but has not yet been split, and =-1 of the node is terminal.
  c A node is terminal if its size is below a threshold value, or if it is
  c all one class, or if all the x-values are equal.  If the current node k
  c is split, then its children are numbered ncur+1 (left), and
  c ncur+2(right), ncur increases to ncur+2 and the next node to be split is
  c numbered k+1.  When no more nodes can be split, buildtree returns to the
  c main program.
 */
/*
      integer a(mdim,nsample),cl(nsample),cat(mdim),
     1     treemap(2,numNodes),bestvar(numNodes),
     1     bestsplit(numNodes), nodestatus(numNodes),ta(nsample),
     1     nodepop(numNodes),nodestart(numNodes),
     1     bestsplitnext(numNodes),idmove(nsample),
     1     ncase(nsample),parent(numNodes),b(mdim,nsample),
     1     jin(nsample),iv(mred),nodeclass(numNodes),mind(mred)
      
      
      double precision tclasspop(nclass),classpop(nclass,numNodes),
     1     tclasscat(nclass,32),win(nsample),wr(nclass),wc(nclass),
     1     wl(nclass),tgini(mdim), xrand
 */          int msplit = 0;
/*      
      zeroInt(nodestatus, numNodes);
      zeroInt(nodestart, numNodes);
      zeroInt(nodepop, numNodes);
      zeroDouble(classpop, nclass * numNodes);
      
      for (i = 0; i < nclass; ++i) classPop[i] = tclassPop[i];
      ncur = 1;
      nodeStart[0] = 1;
      nodePop[0] = *nuse;
      nodeStatus[0] = 2; 
*/
/* 2: not split yet, 1: split, -1: terminal */
      
      /* start main loop */
/*
      for (i = 0; i < numNodes; ++i) {
          if (i > ncur - 1) goto 50;
          if (nodeStatus[kbuild] != 2) goto 30; */
          /* initialize for next call to findbestsplit */
/*          ndstart = nodeStart[kbuild];
          ndend = ndstart + nodePop[kbuild] - 1;
          for (j = 0; j < nclass; ++j) {
              tclasspop[j] = classpop[j, kbuild];
          }
          jstat = 0;
          F77_CALL(findbestsplit)(a,b,cl,mdim,nsample,nclass,cat,ndstart,
     1        ndend,tclasspop,tclasscat,msplit,decsplit,nbest,ncase,
                                  jstat,jin,mtry,win,wr,wc,wl,mred,kbuild,mind);
          if (jstat == 1) { 
              nodeStatus[kbuild] = -1;
              goto 30;
          } else { 
              bestvar[kbuild] = msplit;
              varUsed[msplit - 1] = 1;
              tgini[msplit - 1] += decsplit;
              if (cat[msplit-1] == 1) {
                  bestsplit[kbuild] = a[msplit-1, nbest];
                  bestsplitnext[kbuild] = a[msplit-1, nbest+1];
              } else {
                  bestsplit[kbuild] = nbest;
                  bestsplitnext[kbuild] = 0;
              }
          }
                  
          F77_CALL(movedata)(a,ta,mdim,nsample,ndstart,ndend,idmove,ncase,
                             msplit,cat,nbest,ndendl);
          /-* leftnode no.= ncur+1, rightnode no. = ncur+2. *-/
          nodePop[ncur+1] = ndendl - ndstart + 1;
          nodePop[ncur+2] = ndend - ndendl;
          nodeStart[ncur+1] = ndstart;
          nodeStart[ncur+2] = ndendl + 1;
          /-* find class populations in both nodes *-/
          for (n = ndstart; n <= ndendl; ++n) {
              if (cat[msplit-1] > 1) {
                  nc = ncase[n];
              } else {
                  nc = ncase[n];
              }
              j = class[nc];
              classpop[j,ncur+1] += win[nc];
          }
          for (n=ndendl+1; n <= ndend; ++n) {
              nc = ncase[n];
              j = cl[nc];
              classpop[j,ncur+2] += win[nc];
          }

          /-* check on nodestatus *-/
          nodestatus[ncur + 1] = 2;
          nodestatus[ncur + 2] = 2;
          if (nodepop[ncur + 1] <= ndsize) nodestatus[ncur+1] = -1;
          if (nodepop[ncur + 2] <= ndsize) nodestatus[ncur+2] = -1;
          popt1=0;
          popt2=0;
          for (j=1; j <= nclass; ++j) {
              popt1 += classpop[j,ncur+1];
              popt2 += classpop[j,ncur+2];
          }
          for (j=1; j <= nclass; ++j) {
              if (classpop(j,ncur+1).eq.popt1) nodestatus(ncur+1) = -1;
              if (classpop(j,ncur+2).eq.popt2) nodestatus(ncur+2) = -1;
          }

          treemap(1,kbuild) = ncur + 1;
          treemap(2,kbuild) = ncur + 2;
          parent(ncur+1) = kbuild;
          parent(ncur+2) = kbuild;
          nodestatus(kbuild) = 1;
          ncur += 2;
          if (ncur >= numNodes) break;
      }
      ndbigtree = numNodes;
      for (k=numNodes; k >= 1; --k) {
          if (nodestatus[k] == 0) ndbigtree--;
          if (nodestatus[k] == 2) nodestatus[k] = -1;
      }
      for (k = 1; k <= ndbigtree; ++k) {
          if (nodestatus[k] == -1) {
              pp = 0;
              for (j = 1; j <= nclass; ++j) {
                  if (classpop[j, k] > pp) {
                      nodeclass[k] = j;
                      pp = classpop[j, k];
                  }
                  /-* Break ties at random: *-/
                  if (classpop[j,k] == pp & unif_rand() > 0.5) {
                      nodeclass[k] = j;
                      pp = classpop[j, k];
                  }
              }
          }
      }
}
*/


void predictClassTree(double *x, int n, int mdim, int *doPred, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *cbestsplit, int *bestvar, int *nodeclass,
		      int ndbigtree, int *cat, int nclass,
		      int *jts, int *nodex, int maxcat) {
      
    int icat[32];
    int l, m, ncat, i, j, kt;

    zeroInt(jts, n);
    zeroInt(nodex, n);
    
    /* decode the categorical splits */
    for (i = 0; i < ndbigtree; ++i) {
	if (nodestatus[i] != -1) {
            l = cat[bestvar[i] - 1];
            if (l > 1) {
		ncat = (int) xbestsplit[i];
		unpack(l, ncat, icat);
		for (j = 0; j < l; ++j) {
		    cbestsplit[j + i * maxcat] = icat[j];
		}
            }
	}
    }
    
    for (i = 0; i < n; ++i) {
        if (doPred[i] > 0) continue; /* skip this case */
	kt = 0;
	while (nodestatus[kt] != -1) {
            m = bestvar[kt] - 1;
            if (cat[m] == 1) {
  	        /* Split by a numerical predictor */
	        kt = (x[m + i * mdim] <= xbestsplit[kt]) ?
		    treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    } else {
	        /* Split by a categorical predictor */
	        kt = cbestsplit[(int) x[m + i * mdim] - 1 + kt * maxcat] ?
		    treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    }
	}
	/* Terminal node: assign class label */
	jts[i] = nodeclass[kt];
	nodex[i] = kt + 1;
    }
}
