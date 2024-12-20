// n2pchk.c
#include <stdio.h>
#include "nleq2.h"

void n2pchk(int n, double* x, double* xscal, double* rtol, int* iopt, int* ierr, int liwk, int* iwk, int lrwk, double* rwk) {
    double epmach, small, great;
    zibconst(&epmach, &small);
    great = 1.0 / small;

    // Initialize error code to 0
    *ierr = 0;

    // Check for bad input to dimensional parameter N
    if (n <= 0) {
        *ierr = 20;
        if (iopt[10] >= 1) printf("Bad input to dimensional parameter N.\n");
        return;
    }

    // Check for nonpositive value for RTOL
    if (*rtol <= 0.0) {
        *ierr = 21;
        if (iopt[10] >= 1) printf("Nonpositive value for RTOL supplied.\n");
        return;
    }

    // Check for negative scaling value via vector XSCAL
    for (int i = 0; i < n; i++) {
        if (xscal[i] < 0.0) {
            *ierr = 22;
            if (iopt[10] >= 1) printf("Negative scaling value via vector XSCAL supplied.\n");
            return;
        }
    }

    // Check for invalid fields in IOPT
    for (int i = 0; i < 50; i++) {
        if (iopt[i] < 0) {
            *ierr = 30;
            if (iopt[10] >= 1) printf("One or more fields specified in IOPT are invalid.\n");
            return;
        }
    }

    // Check for sufficient workspace
    int nbroy = (iopt[31] == 1) ? iwk[35] : 0;
    int required_lrwk = (n + nbroy + 15) * n + 61;
    int required_liwk = n + 52;
    if (lrwk < required_lrwk) {
        *ierr = 10;
        if (iopt[10] >= 1) printf("Real workspace too small. Required: %d, Provided: %d\n", required_lrwk, lrwk);
        return;
    }
    if (liwk < required_liwk) {
        *ierr = 10;
        if (iopt[10] >= 1) printf("Integer workspace too small. Required: %d, Provided: %d\n", required_liwk, liwk);
        return;
    }

    // Initialize workspace arrays with zeros
    for (int i = 0; i < liwk; i++) iwk[i] = 0;
    for (int i = 0; i < lrwk; i++) rwk[i] = 0.0;

    // Set the default values for IOPT and workspace arrays
    if (iopt[0] == 0) {
        for (int i = 0; i < 50; i++) iopt[i] = 0;
    }

    // Set machine dependent constants
    iwk[15] = n + 53;
    iwk[16] = (n + 9 + nbroy) * n + 62;

    // Print details if required
    if (iopt[10] >= 2) {
        printf("Internal parameters:\n");
        printf("Starting value for damping factor FCSTART = %.2e\n", rwk[20]);
        printf("Minimum allowed damping factor FCMIN = %.2e\n", rwk[21]);
        printf("Rank-1 updates decision parameter SIGMA = %.2e\n", rwk[22]);
        printf("Initial Jacobian pseudo-rank IRANK = %d\n", iwk[31]);
        printf("Maximum permitted subcondition COND = %.2e\n", rwk[24]);
    }

    // Set lengths of currently required workspaces
    iwk[17] = iwk[15] - 1;
    iwk[18] = iwk[16] - 1;
}