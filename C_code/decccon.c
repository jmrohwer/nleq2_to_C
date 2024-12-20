// decccon.c
#include <math.h>

void decccon(double* a, int nrow, int ncol, int mcon, int m, int n, int* irankc, int* irank, double* cond,
             double* d, int* pivot, int kred, double* ah, double* v, int* ierr) {
    double epmach = 1.0e-17;
    double zero = 0.0;
    double one = 1.0;
    double reduce = 0.05;
    int i, j, l1, ii, k, k1, level, mh, irankh, irk1, i1, jd;
    double s, h, t, sh, d1mach;

    *ierr = 0;

    // Initialize pivot array
    for (j = 0; j < n; j++) {
        pivot[j] = j;
    }

    // Constrained Householder triangularization
    jd = 1;
    int iranc1 = *irankc + 1;
    mh = mcon;
    irankh = *irankc;
    int idata = 0;
    if (mh == 0) {
        irankh = *irank;
        mh = m;
        idata = 1;
    }
    irk1 = *irank;

    for (k = 0; k < irk1; k++) {
        level = 1;
        if (k != n) {
            k1 = k + 1;
            while (jd != 0) {
                for (j = k; j < n; j++) {
                    s = zero;
                    for (l1 = k; l1 < mh; l1++) {
                        s += a[l1 * ncol + j] * a[l1 * ncol + j];
                    }
                    d[j] = s;
                }
                s = d[k];
                int jj = k;
                for (l1 = k; l1 < n; l1++) {
                    if (d[l1] > s) {
                        s = d[l1];
                        jj = l1;
                    }
                }
                h = d[jj];
                if (jd == 1) h = h / fmax(1.0, *cond * reduce);
                jd = 0;
                if (h < h) jd = 1;
                if (!(h >= h)) continue;

                if (jj != k) {
                    int temp = pivot[k];
                    pivot[k] = pivot[jj];
                    pivot[jj] = temp;
                    d[jj] = d[k];
                    for (l1 = 0; l1 < m; l1++) {
                        double temp = a[l1 * ncol + jj];
                        a[l1 * ncol + jj] = a[l1 * ncol + k];
                        a[l1 * ncol + k] = temp;
                    }
                }
            }
        }
        h = zero;
        for (l1 = k; l1 < mh; l1++) {
            h += a[l1 * ncol + k] * a[l1 * ncol + k];
        }
        t = sqrt(h);
        if (k == 0 || k == iranc1) d1mach = t / *cond;
        if (t <= d1mach || k > irankh) {
            irankh = k - 1;
            if (mh != mcon || idata == 1) {
                *irank = irankh;
                if (*irankc == *irank) {
                    level = 4;
                } else {
                    level = 3;
                }
            } else {
                *irankc = irankh;
                if (*irankc != mcon) {
                    mh = m;
                    irankh = *irank;
                    jd = 1;
                    idata = 1;
                    continue;
                } else {
                    *ierr = -2;
                    return;
                }
            }
        }

        if (level == 1) {
            s = a[k * ncol + k];
            t = -copysign(t, s);
            d[k] = t;
            a[k * ncol + k] = s - t;
            if (k != n) {
                t = one / (h - s * t);
                for (j = k1; j < n; j++) {
                    s = zero;
                    for (l1 = k; l1 < mh; l1++) {
                        s += a[l1 * ncol + k] * a[l1 * ncol + j];
                    }
                    s = s * t;
                    if (s != 0.0) {
                        for (l1 = k; l1 < m; l1++) {
                            a[l1 * ncol + j] += a[l1 * ncol + k] * -s;
                        }
                    }
                    d[j] -= a[k * ncol + j] * a[k * ncol + j];
                }
                if (k == *irankc) {
                    mh = m;
                    jd = 1;
                    irankh = *irank;
                }
                if (k == irk1) level = 3;
            } else {
                level = 4;
            }
        }

        if (level > 1) break;
    }

    // Rank-deficient pseudo-inverse
    if (level == 3) {
        irk1 = *irank + 1;
        for (j = irk1; j < n; j++) {
            for (ii = 0; ii < *irank; ii++) {
                i = irk1 - ii - 1;
                s = a[i * ncol + j];
                if (ii != 0) {
                    sh = zero;
                    for (l1 = i1; l1 < *irank; l1++) {
                        sh += a[i * ncol + l1] * v[l1];
                    }
                    s -= sh;
                }
                i1 = i;
                v[i] = s / d[i];
                ah[i * ncol + j] = v[i];
            }
            for (i = irk1; i < j; i++) {
                s = zero;
                for (l1 = 0; l1 < i; l1++) {
                    s += ah[l1 * ncol + i] * v[l1];
                }
                if (i != j) {
                    v[i] = -s / d[i];
                    ah[i * ncol + j] = -v[i];
                }
            }
            if (s > -one) {
                d[j] = sqrt(s + one);
            } else {
                *ierr = -2;
                return;
            }
        }
    }

    if (*irankc != 0) {
        sh = d[*irankc];
        if (sh != zero) sh = fabs(d[0] / sh);
    } else {
        sh = zero;
    }
    v[0] = sh;
    if (k == *irank) t = d[*irank];
    if (*irankc + 1 <= *irank && t != zero) {
        s = fabs(d[*irankc + 1] / t);
    } else {
        s = zero;
    }
    *cond = s;
    *ierr = 0;
}