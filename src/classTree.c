#include <R.h>
#include "rf.h"

void predictClassTree(double *x, int n, int mdim, int *doPred, int *treemap,
		      int *nodestatus, double *xbestsplit,
		      int *cbestsplit, int *bestvar, int *nodeclass,
		      int nrnodes, int ndbigtree, int *cat, int nclass,
		      int *jts, int *nodex, int maxcat) {
      
    int icat[32];
    int k, l, m, ncat, i, j, kt, jcat;

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
	for (k = 0; k < ndbigtree; ++k) {
            if (nodestatus[kt] == -1) {
		/* Terminal node: assign class label */
		jts[i] = nodeclass[kt];
		nodex[i] = kt + 1;
		break;
            }
            m = bestvar[kt] - 1;
            if (cat[m] == 1) {
                /* Split by a numerical predictor */
		kt = (x[m + i * mdim] <= xbestsplit[kt]) ?
		    treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    } else {
		/* Split by a categorical predictor */
		jcat = (int) x[m + i * mdim] - 1;
		kt = cbestsplit[jcat + kt * maxcat] ?
		    treemap[kt * 2] - 1 : treemap[1 + kt * 2] - 1;
	    }
	}
    }
}
