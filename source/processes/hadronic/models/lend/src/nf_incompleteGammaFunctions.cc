/*                            igam.c
 *
 *    Incomplete gamma integral
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                      x
 *                       -
 *                      | |  -t  a-1
 *  igam(a,x)  =        |   e   t   dt.
 *                    | |
 *                     -
 *                      0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30       200000       3.6e-14     2.9e-15
 *    IEEE      0,100      300000       9.9e-14     1.5e-14
 */
/*                            igamc()
 *
 *    Complemented incomplete gamma integral
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                       inf.
 *                         -
 *                        | |  -t  a-1
 *               =        |   e   t   dt.
 *                      | |
 *                       -
 *                        x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 * Tested at random a, x.
 *                a         x                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
 *    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*/

#include "nf_specialFunctions.h"

#if defined __cplusplus
#include <cmath>
#include "G4Exp.hh"
#include "G4Log.hh"
namespace GIDI {
using namespace GIDI;
using namespace std;
#endif

static double big = 4.503599627370496e15;
static double biginv =  2.22044604925031308085e-16;

/*
************************************************************
*/
double nf_incompleteGammaFunctionComplementary( double a, double x, nfu_status *status ) {

    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    *status = nfu_badInput;
    if( !isfinite( x ) ) return( x );
    *status = nfu_Okay;

    if( ( x <= 0 ) || ( a <= 0 ) ) return( 1.0 );
    if( ( x < 1.0 ) || ( x < a ) ) return( nf_gammaFunction( a, status ) - nf_incompleteGammaFunction( a, x, status ) );

    ax = G4Exp( a * G4Log( x ) - x );
    if( ax == 0. ) return( 0.0 );

    if( x < 10000. ) {
        y = 1.0 - a;                        /* continued fraction */
        z = x + y + 1.0;
        c = 0.0;
        pkm2 = 1.0;
        qkm2 = x;
        pkm1 = x + 1.0;
        qkm1 = z * x;
        ans = pkm1 / qkm1;

        do {
            c += 1.0;
            y += 1.0;
            z += 2.0;
            yc = y * c;
            pk = pkm1 * z  -  pkm2 * yc;
            qk = qkm1 * z  -  qkm2 * yc;
            if( qk != 0 ) {
                r = pk / qk;
                t = fabs( ( ans - r ) / r );
                ans = r; }
            else {
                t = 1.0;
            }
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            if( fabs( pk ) > big ) {
                pkm2 *= biginv;
                pkm1 *= biginv;
                qkm2 *= biginv;
                qkm1 *= biginv;
            }
        } while( t > DBL_EPSILON ); } // Loop checking, 11.06.2015, T. Koi
    else {                                  /* Asymptotic expansion. */
        y = 1. / x;
        r = a;
        c = 1.;
        ans = 1.;
        do {
            a -= 1.;
            c *= a * y;
            ans += c;
        } while( fabs( c ) > 100 * ans * DBL_EPSILON ); // Loop checking, 11.06.2015, T. Koi
    }

    return( ans * ax );
}
/*
************************************************************
*/
double nf_incompleteGammaFunction( double a, double x, nfu_status *status ) {
/* left tail of incomplete gamma function:
*
*          inf.      k
*   a  -x   -       x
*  x  e     >   ----------
*           -     -
*          k=0   | (a+k+1)
*/
    double ans, ax, c, r;

    *status = nfu_badInput;
    if( !isfinite( x ) ) return( x );
    *status = nfu_Okay;

    if( ( x <= 0 ) || ( a <= 0 ) ) return( 0.0 );
    if( ( x > 1.0 ) && ( x > a ) ) return( nf_gammaFunction( a, status ) - nf_incompleteGammaFunctionComplementary( a, x, status ) );

    ax = G4Exp( a * G4Log( x ) - x );       /* Compute  x**a * exp(-x) */
    if( ax == 0. ) return( 0.0 );

    r = a;                                                      /* power series */
    c = 1.0;
    ans = 1.0;
    do {
        r += 1.0;
        c *= x / r;
        ans += c;
    } while( c > ans * DBL_EPSILON ); // Loop checking, 11.06.2015, T. Koi

    return( ans * ax / a );
}

#if defined __cplusplus
}
#endif
