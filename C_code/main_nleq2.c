// main_nleq2.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nleq2.h"

// Function prototypes for user-supplied functions
void fcn(int* n, double* x, double* fx, int* iflag);
void jac(int* n, int* ldjac, double* x, double* dfdx, int* iflag);

void zibsec(double* cptim, int* ifail) {
    *cptim = (double)clock() / CLOCKS_PER_SEC;
    *ifail = 0;
}

int main() {
    int irw = 400;
    int iiw = 61;
    int nn = 13;
    int nmaxp = 9;
    int n = 2;
    double eps;
    int iopt[50] = {0};
    int ierr, ifail;
    double x[nn], xscal[nn], rw[irw];
    int iw[iiw];
    double stime, etime, cptime;

    // Opening files
    FILE *data_file = fopen("nleq2.dat", "w");
    FILE *output_file = fopen("nleq2.out", "w");
    if (data_file == NULL || output_file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    printf("Monitor: nleq2.out, Data: nleq2.dat\n");

    while (n <= nmaxp) {
        eps = 1.0e-5;
        int n1 = n + 1;

        for (int i = 0; i < 50; i++) iopt[i] = 0;
        for (int i = 0; i < iiw; i++) iw[i] = 0;
        for (int i = 0; i < irw; i++) rw[i] = 0.0;

        iopt[1] = 1; // Stepwise mode
        iopt[2] = 1; // User supplied Jacobian
        iopt[30] = 3; // Highly nonlinear problem
        iopt[10] = 3; // Print error messages, warnings, and informal messages
        iopt[12] = 9; // Logical unit for error messages
        iopt[13] = 3; // Print iteration monitor
        iopt[14] = 9; // Logical unit for iteration monitor
        iopt[15] = 2; // Print solutions
        iopt[16] = 2; // Logical unit for solutions
        iopt[19] = 1; // Output level for the time monitor
        iopt[20] = 9; // Logical output unit for time monitor

        iw[30] = 200; // Maximum number of iterations

        for (int i = 0; i < n; i++) x[i] = (double)(i + 1) / (double)n1;
        for (int i = 0; i < n; i++) xscal[i] = 0.0;

        ierr = -1;
        int i = 0;

        zibsec(&stime, &ifail);
        while (ierr == -1) {
            nleq2(n, fcn, jac, x, xscal, &eps, iopt, &ierr, iiw, iw, irw, rw);

            // Clear workspace declared not to be used
            int nifree = iw[15];
            for (int k = nifree; k < iiw; k++) iw[k] = 0;
            int nrfree = iw[16];
            for (int k = nrfree; k < irw; k++) rw[k] = 0.0;

            i++;
            fprintf(output_file, "Returned from call %d of NLEQ2\n", i);
        }

        zibsec(&etime, &ifail);
        cptime = etime - stime;
        fprintf(output_file, "\nTime used = %.3f Sec\n", cptime);
        fprintf(output_file, "**********************************************************\n");

        n++;
    }

    fclose(data_file);
    fclose(output_file);
    return 0;
}

// User-supplied function
void fcn(int* n, double* x, double* fx, int* iflag) {
    int i, i1;
    double ti2, ti1, ti, factt;

    for (i = 1; i < *n; i += 2) {
        i1 = i - 1;
        fx[i1] = 0.0;
        fx[i] = (double)(*n) / (double)(i * i - 1);
    }
    if (*n % 2 == 1) {
        fx[*n - 1] = 0.0;
    }
    for (int l = 0; l < *n; l++) {
        factt = 4.0 * x[l] - 2.0;
        ti2 = 1.0;
        ti1 = 0.5 * factt;
        fx[0] = ti1 + fx[0];
        for (int i = 1; i < *n; i++) {
            ti = factt * ti1 - ti2;
            fx[i] = ti + fx[i];
            ti2 = ti1;
            ti1 = ti;
        }
    }
}

// User-supplied Jacobian
void jac(int* n, int* ldjac, double* x, double* dfdx, int* iflag) {
    int i, j;
    double ti2, ti1, ti, factt, tabli2, tabli1, tabli;

    for (j = 0; j < *n; j++) {
        factt = 4.0 * x[j] - 2.0;
        ti2 = 1.0;
        ti1 = 0.5 * factt;
        tabli2 = 0.0;
        tabli1 = 2.0;
        dfdx[j] = tabli1;
        for (i = 1; i < *n; i++) {
            ti = factt * ti1 - ti2;
            tabli = 4.0 * ti1 + factt * tabli1 - tabli2;
            dfdx[i * (*ldjac) + j] = tabli;
            ti2 = ti1;
            ti1 = ti;
            tabli2 = tabli1;
            tabli1 = tabli;
        }
    }
}
