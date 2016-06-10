/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], nf_polevl[];
 *
 * y = nf_polevl( x, coef, N );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as nf_polevl().
 *
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
#include "nf_specialFunctions.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

double nf_polevl( double x, double coef[], int N ) {

    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;

    do {
	    ans = ans * x  +  *p++; }
    while( --i ); // Loop checking, 11.06.2015, T. Koi

    return( ans );
}

/*
************************************************************
*/
/* Evaluate polynomial when coefficient of x^N  is 1.0.  Otherwise same as polevl.  */
double nf_p1evl( double x, double coef[], int N ) {

    double ans;
    double *p;
    int i;

    p = coef;
    ans = x + *p++;
    i = N-1;

    do {
	    ans = ans * x + *p++; }
    while( --i ); // Loop checking, 11.06.2015, T. Koi

    return( ans );
}

#if defined __cplusplus
}
#endif
