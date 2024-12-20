#include <time.h>

void zibsec(double* cptim, int* ifail) {
    // Set CPTIM to cpu time in seconds.
    // This routine is machine dependent.

    *cptim = (double)clock() / CLOCKS_PER_SEC;
    *ifail = 0;
}
