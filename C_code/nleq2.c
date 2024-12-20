// nleq2.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nleq2.h"


void nleq2(int n, void (*fcn)(int*, double*, double*, int*),
           void (*jac)(int*, int*, double*, double*, int*),
           double *x, double *xscal, double *rtol, int *iopt,
           int *ierr, int liwk, int *iwk, int lrwk, double *rwk) {
    // Initialize variables
    int i, j;
    double small, epmach;
    zibconst(&epmach, &small);

    // Workspace initialization
    for (i = 0; i < liwk; i++) iwk[i] = 0;
    for (i = 0; i < lrwk; i++) rwk[i] = 0.0;

    // Default values for IOPT
    if (iopt[0] == 0) {
        for (i = 1; i < 50; i++) iopt[i] = 0;
    }

    // Main NLEQ2 algorithm
    *ierr = 0;
    int mode = iopt[1];
    int jacgen = iopt[2];
    int qsucc = iopt[0];
    int nonlin = iopt[30];

    // Call the parameter checking function
    n2pchk(n, x, xscal, rtol, iopt, ierr, liwk, iwk, lrwk, rwk);
    if (*ierr != 0) return;

    // Scaling
    for (i = 0; i < n; i++) {
        if (xscal[i] < small) xscal[i] = small;
        if (xscal[i] > 1.0 / small) xscal[i] = 1.0 / small;
    }

    // Relative precision
    if (*rtol < epmach * TEN * n) *rtol = epmach * TEN * n;

    int m1 = n;
    int m2 = n;
    int irank = iwk[31];
    double cond = rwk[25];
    int nbroy = iwk[36];

    // Call the core solver N2INT
    n2int(n, fcn, jac, x, xscal, rtol, iwk[31], nonlin, irank, iopt, ierr, lrwk, rwk, iwk[17], n, liwk, iwk, iwk[16], n, m1, m2,
          nbroy, &rwk[51], &rwk[51 + nbroy * n], &rwk[51 + nbroy * n + m1 * n], &rwk[51 + nbroy * n + m1 * n + n], &rwk[51 + nbroy * n + m1 * n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n],
          &rwk[51 + nbroy * n + m1 * n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n + n], mprerr, mprmon, mprsol, luerr, lumon, lusol, niter, ncorr, nfcn, njac, nfcnj, nrejr1, new_, qbdamp);

    // Final scaling
    for (i = 0; i < n; i++) {
        if (xscal[i] < small) xscal[i] = small;
        if (xscal[i] > 1.0 / small) xscal[i] = 1.0 / small;
    }

    // Output final values
    for (i = 0; i < n; i++) {
        x[i] = x[i];
        xscal[i] = xscal[i];
    }

    // Set the achieved accuracy
    *rtol = *rtol;
}