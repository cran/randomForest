#include <R.h>
#include "rf.h"

void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}

void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}

void permuteOOB(int m, double *x, int *in, int nsample, 
		int mdim) {
/* Permute the OOB part of a variable in x. 
 * Argument:
 *   m: the variable to be permuted
 *   x: the data matrix (variables in rows)
 *   in: vector indicating which case is OOB
 *   nsample: number of cases in the data
 *   mdim: number of variables in the data
 */
    double *tp, tmp;
    int i, last, k, nOOB = 0;
    
    tp = (double *) S_alloc(nsample, sizeof(double));

    for (i = 0; i < nsample; ++i) {
	/* make a copy of the OOB part of the data into tp (for permuting) */
	if (in[i] == 0) {
            tp[nOOB] = x[m + i*mdim];
            nOOB++;
        }
    }
    /* Permute tp */
    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
	k = (int) last * unif_rand();
	tmp = tp[last - 1];
	tp[last - 1] = tp[k];
	tp[k] = tmp;
	last--;
    }

    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
	if (in[i] == 0) {
            x[m + i*mdim] = tp[nOOB];
            nOOB++;
	}
    }
}
