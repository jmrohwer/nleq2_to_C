// n2int.c
#include <math.h>
#include "nleq2.h"

// Function prototypes for auxiliary functions
void zibconst(double* epmach, double* small);
void decccon(double* a, int nrow, int ncol, int mcon, int m, int n, int* irankc, int* irank, double* cond,
             double* d, int* pivot, int kred, double* ah, double* v, int* ierr);
void solcon(double* a, int nrow, int ncol, int mcon, int m, int n, double* x, double* b, int irankc, int irank,
            double* d, int* pivot, int kred, double* ah, double* v);

// Numerical Jacobian approximation function
void numerical_jacobian(int n, void (*fcn)(int*, double*, double*, int*), double* x, double* jac, double* fx, double eps) {
    double temp, h;
    int ifail = 0;
    for (int j = 0; j < n; j++) {
        temp = x[j];
        h = eps * fabs(temp);
        if (h == 0.0) h = eps;
        x[j] = temp + h;
        h = x[j] - temp;
        (*fcn)(&n, x, jac + j * n, &ifail);
        x[j] = temp - h;
        h = temp - x[j];
        (*fcn)(&n, x, fx, &ifail);
        for (int i = 0; i < n; i++) {
            jac[i * n + j] = (jac[i * n + j] - fx[i]) / (2.0 * h);
        }
        x[j] = temp;
    }
}

void n2int(int n, void (*fcn)(int*, double*, double*, int*), void (*jac)(int*, int*, double*, double*, int*),
           double *x, double *xscal, double *rtol, int nitmax, int nonlin, int irank, int *iopt, int *ierr,
           int lrwk, double *rwk, int nrfrin, int lrwl, int liwk, int *iwk, int nifrin, int liwl, int m1, int m2,
           int nbroy, double *qa, double *a, double *dx, double *dxq, double *xa, double *xwa, double *fw,
           double *fa, double *eta, double *xw, double *fwk, double *dxqa, double *qu, double *t1, double *t2,
           double *t3, double *fcstrt, double *fcmin, double *sigma, double *sigma2, double *fckeep, double *fc,
           double *fcpri, double *cond, double *dmycor, double *conv, double *sumx, double *dlevf,
           int mprerr, int mprmon, int mprsol, int luerr, int lumon, int lusol, int niter, int ncorr,
           int nfcn, int njac, int nfcnj, int nrejr1, int new_, int qbdamp) {
    int i, j, k;
    int ifail, irankc;
    double epmach, small, great;
    double fcstart, fcbnd, sigma_val, sigma2_val, fcmin_val, cond_val;

    zibconst(&epmach, &small);
    great = 1.0 / small;

    // Initialize variables
    fcstart = *fcstrt;
    fcbnd = 10.0;
    sigma_val = *sigma;
    sigma2_val = *sigma2;
    fcmin_val = *fcmin;
    cond_val = *cond;

    // Initialize counters
    nfcn = 0;
    njac = 0;
    ncorr = 0;
    nrejr1 = 0;
    irankc = irank;

    // Main iteration loop
    for (niter = 0; niter < nitmax; niter++) {
        // Function evaluation
        (*fcn)(&n, x, fw, &ifail);
        nfcn++;
        if (ifail < 0) {
            *ierr = 82;
            iwk[22] = ifail;
            return;
        }

        // Check convergence
        *conv = 0.0;
        for (i = 0; i < n; i++) {
            *conv += (fw[i] / xscal[i]) * (fw[i] / xscal[i]);
        }
        *conv = sqrt(*conv / (double)n);
        if (*conv <= *rtol) {
            *ierr = 0;
            return;
        }

        // Jacobian evaluation
        if (iopt[2] == 1) {
            (*jac)(&n, &n, x, a, &ifail);
            njac++;
            if (ifail < 0) {
                *ierr = 83;
                iwk[22] = ifail;
                return;
            }
        } else {
            // Numerical Jacobian approximation
            numerical_jacobian(n, fcn, x, a, fw, epmach);
            njac++;
        }

        // Solve the linear system
        decccon(a, n, n, 0, n, n, &irankc, &irank, &cond_val, rwk, iwk, 0, qa, xwa, &ifail);
        if (ifail != 0) {
            *ierr = 80;
            iwk[22] = ifail;
            return;
        }
        solcon(a, n, n, 0, n, n, dx, fw, irankc, irank, rwk, iwk, 0, qa, xwa);

        // Line search and update
        for (k = 0; k < 10; k++) {
            for (i = 0; i < n; i++) {
                x[i] = xa[i] - fcstart * dx[i];
            }

            (*fcn)(&n, x, fa, &ifail);
            nfcn++;
            if (ifail < 0) {
                *ierr = 82;
                iwk[22] = ifail;
                return;
            }

            *sumx = 0.0;
            for (i = 0; i < n; i++) {
                *sumx += (fa[i] / xscal[i]) * (fa[i] / xscal[i]);
            }
            *sumx = sqrt(*sumx / (double)n);
            if (*sumx <= fcbnd * (*conv)) {
                break;
            }

            fcstart *= 0.5;
        }
    }

    // If max iterations reached without convergence
    *ierr = 2;
}
