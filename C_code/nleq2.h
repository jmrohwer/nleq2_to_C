// nleq2.h
#ifndef NLEQ2_H
#define NLEQ2_H

#ifdef __cplusplus
extern "C" {
#endif

void nleq2(int n, void (*fcn)(int*, double*, double*, int*),
           void (*jac)(int*, int*, double*, double*, int*),
           double *x, double *xscal, double *rtol, int *iopt,
           int *ierr, int liwk, int *iwk, int lrwk, double *rwk);

#ifdef __cplusplus
}
#endif

#endif // NLEQ2_H