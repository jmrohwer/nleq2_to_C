// solcon.c
#include <math.h>

void solcon(double* a, int nrow, int ncol, int mcon, int m, int n, double* x, double* b, int irankc, int irank,
            double* d, int* pivot, int kred, double* ah, double* v) {
    double zero = 0.0;
    int i, ii, l1, l2, j, j1, i1;
    double s, sh;

    if (irank == 0) {
        for (l1 = 0; l1 < n; l1++) {
            x[l1] = zero;
        }
        return;
    }
    if (irank <= irankc && irank != n) {
        int iranc1 = irankc + 1;
        for (l1 = iranc1; l1 < n; l1++) {
            v[l1] = zero;
        }
    }
    if (kred >= 0 && (m != 1 || n != 1)) {
        int mh = mcon;
        if (irankc == 0) mh = m;
        for (int j = 0; j < irank; j++) {
            s = zero;
            for (l1 = j; l1 < mh; l1++) {
                s += a[l1 * ncol + j] * b[l1];
            }
            s = s / (d[j] * a[j * ncol + j]);
            for (l1 = j; l1 < m; l1++) {
                b[l1] += a[l1 * ncol + j] * s;
            }
            if (j == irankc) mh = m;
        }
    }

    int irk1 = irank + 1;
    for (ii = 0; ii < irank; ii++) {
        i = irk1 - ii - 1;
        i1 = i + 1;
        s = b[i];
        if (ii != 0) {
            sh = zero;
            for (l1 = i1; l1 < irank; l1++) {
                sh += a[i * ncol + l1] * v[l1];
            }
            s -= sh;
        }
        v[i] = s / d[i];
    }
    if (irank != n && irank != irankc) {
        for (int j = irk1; j < n; j++) {
            s = zero;
            for (l1 = 0; l1 < j; l1++) {
                s += ah[l1 * ncol + j] * v[l1];
            }
            v[j] = -s / d[j];
        }
        for (int jj = 0; jj < n; jj++) {
            j = n - jj - 1;
            s = zero;
            if (jj != 0) {
                for (l1 = j1; l1 < n; l1++) {
                    s += ah[j * ncol + l1] * v[l1];
                }
            }
            if (jj != 0 && j <= irank) {
                v[j] = v[j] - s;
            } else {
                j1 = j;
                v[j] = -(s + v[j]) / d[j];
            }
        }
    }

    for (l1 = 0; l1 < n; l1++) {
        l2 = pivot[l1];
        x[l2] = v[l1];
    }
}
