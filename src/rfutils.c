#include <R.h>
#include "rf.h"

void zeroInt(int *x, int length) {
	int i;
	for (i=0; i < length; ++i) x[i] = 0;
	/* memset(x, 0, length * sizeof(int)); */
}

void zeroDouble(double *x, int length) {
	int i;
	for (i=0; i < length; ++i) x[i] = 0.0;
	/* memset(x, 0, length * sizeof(double)); */
}

void createClass(double *x, int realN, int totalN, int mdim) {
/* Create the second class by bootstrapping each variable independently. */
    int i, j, k;
    for (i = realN; i < totalN; ++i) {
        for (j = 0; j < mdim; ++j) {
            k = (int) unif_rand() * realN;
            x[j + i * mdim] = x[j + k * mdim];
        }
    }
}

void prepare(int *cl, const int nsample, const int nclass, const int ipi, 
	     double *pi, double *pid, int *nc, double *wtt) {
    int i, j;
    double sump=0.0;
    
    zeroInt(nc, nclass);
    for (i = 0; i < nsample; ++i) nc[cl[i] - 1] ++;

    if (ipi == 0) {
       for (j = 0; j < nclass; ++j) pi[j] = ((double) nc[j]) / nsample;
    }   
    for (i = 0; i < nclass; ++i) sump += pi[i];
    for (i = 0; i < nclass; ++i) pi[i] = pi[i] / sump;

    for (i = 0; i < nclass; ++i) {
        pid[i] = nc[i] ? pi[i] * nsample / nc[i] : 0.0;
    }
    for (j = 0; j < nsample; ++j) wtt[j] = pid[cl[i]-1];
}

void makeA(double *x, const int mdim, const int nsample, int *cat, int *a, 
           int *b) {
    /* makeA() constructs the mdim by nsample integer array a.  if there are
     * less than 32,000 cases, this can be declared integer*2, otherwise
     * integer*4. For each numerical variable with values x(m,n),
     * n=1,...,nsample, the x-values are sorted from lowest to highest.
     * Denote these by xs(m,n).  Then a(m,n) is the case number in which
     * xs(m,n) occurs. The b matrix is also contructed here.
     * if the mth variable is categorical, then a(m,n) is the category of the
     * nth case number.
     */
    int i, j, n1, n2, *index;
    double *v;

    v     = (double *) Calloc(nsample, double);
    index = (int *) Calloc(nsample, int);

    for (i = 0; i < mdim; ++i) {
        if (cat[i] == 1) { /* numerical predictor */
            for (j = 0; j < nsample; ++j) {
                v[j] = x[i + j * mdim];
                index[j] = j + 1;
            }
            R_qsort_I(v, index, 1, nsample);

            /*  this sorts the v(n) in ascending order. isort(n) is the case 
                number of that v(n) nth from the lowest (assume the original 
                case numbers are 1,2,...).  */
            for (j = 0; j < nsample-1; ++j) {
                n1 = index[j];
                n2 = index[j + 1];
                a[i + j * mdim] = n1;
                if (j == 0) b[i + (n1-1) * mdim] = 1;
                b[i + (n2-1) * mdim] =  (v[j] < v[j + 1]) ?
                    b[i + (n1-1) * mdim] + 1 : b[i + (n1-1) * mdim];
            }
            a[i + (nsample-1) * mdim] = index[nsample-1];
        } else { /* categorical predictor */
            for (j = 0; j < nsample; ++j) 
                a[i + j*mdim] = (int) x[i + j * mdim];
        }
    }
    Free(index);
    Free(v);
}


void modA(int *a, int *nuse, const int nsample, const int mdim,
	  int *cat, const int maxcat, int *ncase, int *jin) {
    int i, j, k, m, nt;

    *nuse = 0;
    for (i = 0; i < nsample; ++i) if (jin[i]) (*nuse)++;
    
    for (i = 0; i < mdim; ++i) {
      k = 0;
      nt = 0;
      if (cat[i] == 1) {
          for (j = 0; j < nsample; ++j) {
              if (jin[a[i + k * mdim] - 1]) {
                  a[i + nt * mdim] = a[i + k * mdim];
                  k++;
              } else {
                  for (m = 0; m < nsample - k; ++m) {
                      if (jin[a[i + (k + m) * mdim] - 1]) {
                          a[i + nt * mdim] = a[i + (k + m) * mdim];
                          k += m + 1;
                          break;
                      }
                  }
              }
              nt++;
              if (nt >= *nuse) break;
          }
      }
    }
    if (maxcat > 1) {
        k = 0;
        nt = 0;
        for (i = 0; i < nsample; ++i) {
            if (jin[k]) {
                k++;
                ncase[nt] = k;
            } else {
                for (j = 0; j < nsample - k; ++j) {
                    if (jin[k + j]) {
                        ncase[nt] = k + j + 1;
                        k += j + 1;
                        break;
                    }
                }
            }
            nt++;
            if (nt >= *nuse) break;
        }
    }
}

void Xtranslate(double *x, int mdim, int nrnodes, int nsample, 
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize) {
/*
  c this subroutine takes the splits on numerical variables and translates them
  c back into x-values.  It also unpacks each categorical split into a 
  c 32-dimensional vector with components of zero or one--a one indicates that 
  c the corresponding category goes left in the split.
*/

    int i, m;

    for (i = 0; i < treeSize; ++i) {
	if (nodestatus[i] == 1) {
	    m = bestvar[i] - 1;
	    if (cat[m] == 1) {
		xbestsplit[i] = 0.5 * (x[m + (bestsplit[i] - 1) * mdim] +
				       x[m + (bestsplitnext[i] - 1) * mdim]);
	    } else {
		xbestsplit[i] = (double) bestsplit[i];
	    }
	}
    }
}

void permuteOOB(int m, double *x, int *in, int nsample, int mdim) {
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
    
    tp = (double *) Calloc(nsample, double);

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

    /* Copy the permuted OOB data back into x. */
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
	if (in[i] == 0) {
            x[m + i*mdim] = tp[nOOB];
            nOOB++;
	}
    }
    Free(tp);
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
